
/*
 *  segemehl.c
 *  segemehl
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 02:50:57 PM CEST
 *
 *  Revision of last commit: 
 *  $Rev: 103 $
 *  $Author: steve $
 *  $Date: 2008-12-10 15:18:18 +0100 (Wed, 10 Dec 2008) $
 *
 *
 *  $Id: segemehl.c 103 2008-12-10 14:18:18Z steve $
 *  $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "memory.h"
#include "biofiles.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include <sys/types.h>
#include <unistd.h>
#include <sys/times.h>
#include "vtprogressbar.h"
#include "manout.h"
#include <time.h>
#include "sort.h"
#include "list.h"
#include "biofiles.h"
#include "kdmatch.h"
#include "debug.h"
#include "info.h"
#include <pthread.h>
#include "citation.h"
#include "kdseed.h"
#include "manopt.h"
#include "segemehl.h"

/*//global variables
 unsigned int  *EDtabcolumn,   // a column of the Edit Distance table
       *Rtabcolumn,    // a column of the Row index table
       *Ctab;    // a table to crosspoints
 
 unsigned int  rekCounter=0;
 int  **EMatrix;
 int  *LinEMatrix;
*/
pthread_mutex_t updatemtx;
unsigned char mute=0;

void*
checkclock(void *args) {
  checkthreadinfo_t *t;

  sleep(2);
  cursorVisible();
  t = (checkthreadinfo_t*) args;
  initProgressBarVT();

  while (pthread_mutex_trylock(&updatemtx) != 0) {
    progressBarVT("reads matched.", t->noofseqs, (*t->counter), 25);
  }

  cursorVisible();
  fprintf(stderr, "\n");
  return NULL;
}


void*
kmismatchworker(void *args) {

  segemehlinfo_t *t;
  t = (segemehlinfo_t*) args;
  kmismatch(t->space, t->sarray, t->reads, t->kMis, t->counter, 
      t->rep_type, 0, t->dev);

  return NULL;
}

void*
kdseedworker(void *args) {
  segemehlinfo_t *t;
  t = (segemehlinfo_t*) args;
  
  matchkdseed(t->space, t->sarray, t->reads, 
      t->minsize,  t->outfile, t->counter, 0, t->s_ext, t->p_mis, t->Xoff, 
      t->k_p, t->rep_type, t->bestonly, t->bedist, t->align, t->maxevalue, t->accuracy, t->M, 
      t->matchingstat, t->dev, t->uninformativedev);

  return NULL;
}


int main(int argc, char** argv) {
  
  segemehlinfo_t info, *th_info;
  manopt_arg *unflagged;
  manopt_optionset optset;
  manopt_intconstraint threadconstraint;
  manopt_intconstraint accuracyconstraint;
  
  int *space = NULL;
  int i =0;
  unsigned int *suflinktable,
               counter=0; 
  unsigned char index = 0,
                kmismatchmode =0,
                brief=0;
  double difsuf,
         difmatch;
  time_t startsuf, endsuf;
  time_t startmatch, endmatch;
  time_t rawtime;
  struct tm *timeinfo;

  fasta_t  **chopsuey;  
  pthread_t *threads;
  pthread_t clockthread;
  checkthreadinfo_t ch_info;
  manopt_arg *dbfilenames;
  
  threadconstraint.max = 3000;
  threadconstraint.min = 1;
  accuracyconstraint.max=100;
  accuracyconstraint.min=0;
 
  setdefault(&info);
  manopt_initoptionset(&optset, argv[0], NULL, 
      "  Heuristic mapping of short sequences\n",
      "  SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      VERSION ,
      "  Please report bugs to steve@bioinf.uni-leipzig.de");

  manopt_blockseparator(&optset, "INPUT");
  manopt(&optset, LISTOPT, 1, 'd', "database", 
      "list of path/filename(s) of database sequence(s)", "<file> [<file>]", 
       NULL, NULL);
  manopt(&optset, REQSTRINGOPT, 0, 'q', "query", 
      "path/filename of query sequence", "<file>", NULL, &info.queryfilename);
  manopt(&optset, REQSTRINGOPT, 0, 'i', "index", 
      "path/filename of db index", "<file>", NULL, &info.idxfilename);
  manopt(&optset, REQSTRINGOPT, 0, 'x', "generate", 
      "generate db index and store to disk", "<file>", NULL, &info.idxfilename);
  manopt_blockseparator(&optset, "GENERAL");
  manopt(&optset, REQUINTOPT, 0, 'm', "minsize", 
      "minimum size of queries", "<n>", NULL, &info.minsize);
  manopt(&optset, FLAG, 0, 's', "silent", 
      "shut up!", NULL, NULL, &mute);
  manopt(&optset, FLAG, 0, 'b', "brief", 
      "brief output", NULL, NULL, &brief);
  manopt(&optset, FLAG, 0, 'c', "checkidx", 
      "check index", NULL, NULL, &info.check);
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
      "start <n> threads", "<n>", &threadconstraint, &info.threadno);
  manopt(&optset, REQSTRINGOPT, 0, 'o', "outfile", 
      "dropoff parameter for extension", "<n>", NULL, &info.outfile);
  manopt(&optset, REQSTRINGOPT, 0, 'u', "nomatchfilename", 
      "filename for unmatched reads", "<file>", NULL, &info.uninformativename); 
  manopt_blockseparator(&optset, "SEEDPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'D', "differences", 
      "search seeds initially with <n> differences", "<n>", NULL, &info.k_p);
  manopt(&optset, REQDBLOPT, 0, 'E', "evalue", 
      "max evalue", "<double>", NULL, &info.maxevalue);  
  manopt(&optset, REQUINTOPT, 0, 'M', "maxinterval", 
      "maximum width of a suffix array interval, i.e. a query seed will be omitted if it matches more than <n> times", "<n>", NULL, &info.M); 
  manopt_blockseparator(&optset, "SEEDEXTENSIONPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'e', "extensionscore", 
      "score of a match during extension", "<n>", NULL, &info.s_ext); 
  manopt(&optset, REQUINTOPT, 0, 'p', "extensionpenalty", 
      "penalty for a mismatch during extension", "<n>", NULL, &info.p_mis);
  manopt(&optset, REQUINTOPT, 0, 'X', "dropoff", 
      "dropoff parameter for extension", "<n>", NULL, &info.Xoff);
  manopt_blockseparator(&optset, "ALIGNPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'A', "accuracy", 
      "min percentage of matches per read in semi-global alignment", "<n>", NULL, &info.accuracy);
  manopt(&optset, REQUINTOPT, 0, 'H', "hitstrategy", 
      "report only best scoring hits (=2), best hits per strand (=1) or all (=0)", NULL, NULL, &info.bestonly);  
  manopt(&optset, FLAG, 0, 0, "showalign", 
      "show alignments", NULL, NULL, &info.align);
/*  manopt_blockseparator(&optset, "ALTERNATIVE MODES");
  manopt(&optset, REQUINTOPT, 0, 'k', "kmismatch", 
      "search queries allowing <n> mismatches", "<n>", NULL, &info.kMis);
  manopt(&optset, FLAG, 0, 'S', "matchstat", 
      "calculate matching statistics for queries", "<n>", NULL, 
      &info.matchingstat); 
  manopt(&optset, REQUINTOPT, 0, 'r', "report", 
      "select different report type", "<n>", NULL, &info.rep_type); 
  manopt(&optset, REQUINTOPT, 0, 'B', "bestedist", 
      "report hits not more than <n> errors from best edist away", "<n>", NULL, &info.bedist);
  manopt(&optset, REQUINTOPT, 0, 'F', "fusion", 
      "report potential fusion transcripts if at least <n> bp map to different sites", "<n>", NULL, &info.fusion);
*/  
  unflagged = manopt_getopts(&optset, argc, argv);
  index = manopt_isset(&optset, 'x', NULL);

  if(!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
    manopt_help(&optset, "please give index filename using -i XOR -x option\n");
  } else if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unkown argument(s)\n");
  }
  if(manopt_isset(&optset, 'k', "kmismatch")) {
    kmismatchmode = 1;
  }
  if(manopt_isset(&optset, 'b', "brief")) {
    info.rep_type = 5;
  }

  if(info.outfile) {
    info.dev=fopen(info.outfile, "w");
  }

  if(info.uninformativename != NULL)
  info.uninformativedev = fopen(info.uninformativename, "w");
  
  if (info.dev == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", info.outfile);
    exit(EXIT_FAILURE);
  }

  dbfilenames = manopt_getarg(&optset, 'd', "database");
  info.fasta = getfastaset(space, dbfilenames->values, 
      dbfilenames->noofvalues, 1, 0);

  for(i=0; i < info.fasta->noofseqs; i++) {
    info.totallength += info.fasta->seqs[i]->length; 
  }

  if(index) {
    time (&startsuf);
    info.sarray = constructSufArr(space, info.fasta->seqs, 
        info.fasta->noofseqs, NULL, mute); 

  if (info.check) {
    DBG("checking suffixarray %s\n", info.idxfilename);
    for(i=1; i < info.sarray->numofsuffixes-1; i++) {
      if(!mute) {
        progressBarVT("suffixes checked.", info.sarray->numofsuffixes, i, 25);
      }
      if (strcmp(&info.sarray->seq->sequences[info.sarray->suftab[i-1]],
            &info.sarray->seq->sequences[info.sarray->suftab[i]]) > 0) {
        DBG("suffix array '%s' corrupted! Exit forced.\n", info.idxfilename);
      }
    }
  }
    MSG("constructing lcp.\n");
    constructLcp(space, info.sarray);
    MSG("deleting inv_suftab\n");
    FREEMEMORY(space, info.sarray->inv_suftab);
    info.sarray->inv_suftab = NULL;

    MSG("constructing child tab.\n");
    constructchildtab(space, info.sarray);
    MSG("constructing suffix links.\n");
    MSG("constructing id.\n");
    computeId(space, info.sarray);
    MSG("constructing suflinks - bottom up.\n");
    suflinktable = getsufsucc(space, info.sarray);
    MSG("constructing suflinks - top down.\n");
    constructsuflinks(space, info.sarray, suflinktable);
    FREEMEMORY(space, suflinktable);
    time (&endsuf);
    difsuf = difftime (endsuf, startsuf);
    NFO("building  the suffix array has taken %f seconds.\n", difsuf);
    NFO("total length of suffix array was %u.\n", info.totallength);

  } else {

    time (&startsuf);
    NFO("reading suffix array '%s' from disk.\n", info.idxfilename);
    info.sarray=readSuffixarray(space, info.idxfilename, info.fasta->seqs, 
        info.fasta->noofseqs, mute); 
    time (&endsuf);
    difsuf = difftime (endsuf, startsuf);

    NFO("reading the suffix array has taken %f seconds.\n", difsuf);
  }

  if (info.check) {
    DBG("checking suffixarray %s\n", info.idxfilename);
    for(i=1; i < info.sarray->numofsuffixes-1; i++) {
      if(!mute) {
        progressBarVT("suffixes checked.", info.sarray->numofsuffixes, i, 25);
      }
      if (strcmp(&info.sarray->seq->sequences[info.sarray->suftab[i-1]],
            &info.sarray->seq->sequences[info.sarray->suftab[i]]) > 0) {
        DBG("suffix array '%s' corrupted! Exit forced.\n", info.idxfilename);
      }
    }
  }

  if(index && info.idxfilename) {
    NFO("writing suffix '%s' array to disk.\n", info.idxfilename);
    writeSuffixarray(info.sarray, info.idxfilename); 
  }

  if(info.queryfilename) {
    info.reads = readfasta(space, info.queryfilename, 1, 0);
    time(&rawtime);
    timeinfo = localtime (&rawtime);

    fprintf(info.dev, "#segemehl\t%s", asctime(timeinfo));
    fprintf(info.dev, "#query: %s\n", info.queryfilename);
    fprintf(info.dev, "#subject: %s\n",info.idxfilename);

    if(info.kMis == 0) {
      fprintf(info.dev, "#minsize=%d, diff_seed=%d, score_ext=%d, ", 
          info.minsize, info.k_p, info.s_ext);
      fprintf(info.dev, "penalty_mis:%d, Xoff=%d, acc=%d, maxEvalue: %.5f, hitstrategy: %d \n", 
          info.p_mis, info.Xoff, info.accuracy, info.maxevalue, info.bestonly);
      matchHeader(info.dev, info.rep_type); 
    } else {
      fprintf(info.dev, "#k-mismatch mode with k=%d\n", info.kMis);
    }

    if (info.threadno > 1) {
      if (info.threadno > info.reads->noofseqs || info.reads->noofseqs < 2) {
        DBG("%d threads for <100 seqs req. Exit forced\n", info.threadno);
        exit(EXIT_FAILURE);
      }
      info.counter=&counter;
      NFO("starting %d threads.\n", info.threadno);
      chopsuey = chopfasta(space, info.reads, info.threadno);
      th_info = ALLOCMEMORY(space, NULL, segemehlinfo_t, info.threadno);
      threads = ALLOCMEMORY(space, NULL, pthread_t, info.threadno);
      ch_info.noofseqs = info.reads->noofseqs;
      ch_info.counter = &counter;

      if (!mute) {
        pthread_mutex_lock(&updatemtx);
        pthread_create(&clockthread, NULL, checkclock, &ch_info);
      }

      for(i=0; i < info.threadno; i++) {
        NFO("%d reads in thread %d.\n", chopsuey[i]->noofseqs, i);
      }

      time (&startmatch);

      for(i=0; i < info.threadno; i++) {
        memmove(&th_info[i], &info, sizeof(segemehlinfo_t));
        th_info[i].reads = chopsuey[i];
        if (!kmismatchmode) {
          pthread_create(&threads[i], NULL, kdseedworker, &th_info[i]);
        } else {
          pthread_create(&threads[i], NULL, kmismatchworker, &th_info[i]);
        }
      }

      for(i=0; i < info.threadno; i++) {
        pthread_join(threads[i], NULL); 
      } 

      if(!mute) {
        /*notification via mutex - why use signals?*/
        pthread_mutex_unlock(&updatemtx);
        pthread_join(clockthread, NULL);
      }
      fflush(info.dev);
      time (&endmatch);
      difmatch = difftime (endmatch, startmatch);
      NFO("threaded matching w/ suffixarray has taken %f seconds.\n", difmatch);
      
      for (i=0; i < info.threadno; i++) {
        free(chopsuey[i]->seqs);
        free(chopsuey[i]);
      }
      FREEMEMORY(space, chopsuey);
      FREEMEMORY(space, th_info);
      FREEMEMORY(space, threads);
    } else {
      if (!kmismatchmode) { 
        time (&startmatch);
        matchkdseed(space, info.sarray, info.reads, info.minsize, 
            info.outfile, 0, mute, info.s_ext, info.p_mis, info.Xoff, 
            info.k_p, info.rep_type, info.bestonly, info.bedist, info.align, info.maxevalue, 
            info.accuracy, info.M, info.matchingstat, info.dev, info.uninformativedev);
        time (&endmatch);
        difmatch = difftime (endmatch, startmatch);
        NFO("matching w/ suffixarray has taken %f seconds.\n", difmatch); 
      } else {
        time (&startmatch);
        kmismatch(info.space, info.sarray, info.reads, info.kMis, 
            NULL, info.rep_type, mute, info.dev);
        time (&endmatch);
        difmatch = difftime (endmatch, startmatch);
        NFO("matching w/ suffixarray has taken %f seconds.\n", difmatch); 
      }  
    }

    destructfasta(space, info.reads);
    FREEMEMORY(space, info.reads);
  }

  if (info.outfile) {
    fclose(info.dev);
  }
  
  if(info.uninformativename != NULL)
  fclose(info.uninformativedev);

  destructfasta(space, info.fasta);
  destructSufArr(space, info.sarray);
  FREEMEMORY(space, info.fasta);

  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  free(unflagged);

  NFO("\nGoodbye.\n %s\n", citerand());  
  return EXIT_SUCCESS;

}

