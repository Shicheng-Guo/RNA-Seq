
/*
 *  matfile.c
 *  match files
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/25/2010 03:42:37 PM CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "bitVector.h"
#include "matchfiles.h"
#include "browsematchfiles.h"
#include "matfile.h"
#include "iupac.h"
#include "info.h"
#include "manopt.h"
#include "evalmatchfiles.h"
#include "splicesites.h"
#include "startsites.h"

unsigned char mute = 0;
char *ntcode;


/*---------------------------------- test -----------------------------------
 *    
 * @brief get matchfile stats
 * @author Steve Hoffmann 
 *   
 */
 
  void
teststats (void *space, matchfile_t *file, fasta_t *fasta)
{

  //double ws[]={1, 1};
  Uint ret;
  matchfileSampleStats_t *stats;

  stats = bl_matchfileInitSampleStats(space, 1000000, 100, 6, 10, .15);
  //bl_matchfileDumpFileStats(space, file);
  //bl_matchfileEvalCrossSections(space, file, fasta, NULL, NULL);
  //file->index->stats = stats;

  MSG("evaluating error density.\n");
  //ret = bl_matchfileSampleCrossSections(space, file, fasta, 10000, 
  //    bl_matchfileGetErrorDensity, stats);

  //if(!ret) return;

  //stats->e_ll=gmm(space, stats->e, stats->e_N, 1, stats->e_mu, 
  //    stats->e_sd, ws, 2, 10000);

  //stats->px =stats->e_mu[0] + stats->e_sd[0];
  //while(univarnormcdf(stats->px, stats->e_mu[0], stats->e_sd[0]) < 0.99) 
  //  stats->px += 0.00001;

  MSG("sample cross sections.\n");
  ret = bl_matchfileSampleCrossSections(space, file, fasta, 100000, 
      bl_matchfileGetConditionals, stats);

  if(!ret) return;

  //stats->V_ll=gmm(space, stats->V, stats->V_N, 1, 
  //    &stats->V_mu, &stats->V_sd, ws, 1, 100000); 

  //stats->Vx_ll=gmm(space, stats->Vx, stats->Vx_N, 1, 
  //    &stats->Vx_mu, &stats->Vx_sd, ws, 1, 100000);


  return ;
}



/*---------------------------------- stats -----------------------------------
 *    
 * @brief get matchfile stats
 * @author Steve Hoffmann 
 *   
 */
 
  void
getstats (void *space, matchfile_t *file, fasta_t *fasta)
{

  double ws[]={1, 1};
  Uint ret;
  matchfileSampleStats_t *stats;

  MSG("evaluating statistics\n");
  stats = bl_matchfileInitSampleStats(space, 1000000, 100, 6, 10, .15);
  //bl_matchfileDumpFileStats(space, file);
  //bl_matchfileEvalCrossSections(space, file, fasta, NULL, NULL);
  file->index->stats = stats;

  MSG("evaluating error density.\n");
  ret = bl_matchfileSampleCrossSections(space, file, fasta, 10000, 
      bl_matchfileGetErrorDensity, stats);

  if(!ret) return;

  stats->e_ll=gmm(space, stats->e, stats->e_N, 1, stats->e_mu, 
      stats->e_sd, ws, 2, 10000);

  stats->px =stats->e_mu[0] + stats->e_sd[0];
  while(univarnormcdf(stats->px, stats->e_mu[0], stats->e_sd[0]) < 0.99) 
    stats->px += 0.00001;

  MSG("sample cross sections.\n");
  ret = bl_matchfileSampleCrossSections(space, file, fasta, 100000, 
      bl_matchfileGetConditionals, stats);

  if(!ret) return;

  stats->V_ll=gmm(space, stats->V, stats->V_N, 1, 
      &stats->V_mu, &stats->V_sd, ws, 1, 100000); 

  stats->Vx_ll=gmm(space, stats->Vx, stats->Vx_N, 1, 
      &stats->Vx_mu, &stats->Vx_sd, ws, 1, 100000);


  return ;
}

/*----------------------------------- view -----------------------------------
 *    
 * @brief view the matchfiles
 * @author Steve Hoffmann 
 *   
 */

void
view (void *space, matchfile_t **files, Uint nooffiles, fasta_t *fasta, 
    annotationtrack_t *bed)
{
 
  MSG("starting viewer.\n");
  bl_matchfileViewer(space, files, nooffiles, fasta, bed, 100, 1000);
  return ;
}

/*----------------------------------- eval -----------------------------------
 *    
 * @brief evaluation of matchfiles
 * @author Steve Hoffmann 
 *   
 */
 
void
eval (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta)
{
 
  MSG("evaluating x-sections.\n");
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, 
      bl_matchfileTest, NULL);
 
  return ;
}

/*-------------------------------- evalstarts --------------------------------
 *    
 * @brief get start sites
 * @author Steve Hoffmann 
 *   
 */
 
void
evalstarts (void *space, matchfile_t **files, int *groups, Uint nooffiles, 
    fasta_t *fasta, annotationtrack_t *bed )
{
 
  Uint *cntr;
//  Uint i;
  cntr = ALLOCMEMORY(space, NULL, Uint, 255);
  memset(cntr, 0, sizeof(Uint)*255);
 
  MSG("start sites.\n");
/*
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, 
      bl_matchfileStartSites, cntr);

  for(i=0; i < 255; i++) {
    printf("%d\t%d\n", i, cntr[i]);
  }
*/
  
  MSG("coverage.\n");
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, 
      bl_coverage, cntr);

  return ;
}

/*-------------------------------- writewiggle -------------------------------
 *    
 * @brief get start sites
 * @author Steve Hoffmann 
 *   
 */
 
void
writeexpression (void *space, matchfile_t **files, int *groups, Uint nooffiles, 
    fasta_t *fasta, annotationtrack_t *bed)
{
 
  Uint i;
  char *filename;
  MSG("write expression file.\n");

  for(i=0; i < nooffiles; i++) { 
    filename = ALLOCMEMORY(space, NULL, char, strlen(files[i]->filename)+5);
    memmove(filename, files[i]->filename, strlen(files[i]->filename));
    sprintf(&filename[strlen(files[i]->filename)], ".wig");

    bl_writeexpression(space, filename, strlen(filename), files[i], fasta, 
        10000000, 1000000);
  }

  return ;
}



/*---------------------------------- splice ----------------------------------
 *    
 * @brief dump the splice sites
 * @author Steve Hoffmann 
 *   
 */
 
void
evalsplice (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, 
    annotationtrack_t *bed, char *filename, char *transfilename, Uint minsplitno)
{
  Uint i, cidx, pos, len;
  char *chr, *chr2, *ptr, *base;
  splitmap_t map;
  splicemap_t *sm;
  FILE *fp1, *fp2;

  bl_matchfileInitSplitMap(space, &map, bed, files, fasta);
  MSG("eval splice sites.\n");
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, 
      bl_matchfileSplit, &map);
  MSG("condensing sites.\n");
  sm = bl_matchfileSpliceMap (space, &map, 11, minsplitno);
  MSG("writing splice sites to stdout.\n");



  fp1 = fopen(filename, "w");
  if (fp1 == NULL) {
    fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", filename);
    exit(-1);
  } 


  fp2 = fopen(transfilename, "w");
  if (fp2 == NULL) {
    fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", transfilename);
    exit(-1);
  } 


  base = bl_basename(filename);

  //printsplice(space, sm, stdout);
  printsplicebed(space, sm, base, fp1, fp2);
 
  MSG("writing splice sites to gff.\n");
  if(bed) {
    bl_annotationtrackGetStats (space, bed);
    bl_matchfileSpliceAnnotation(space, sm, bed);
    bl_GFFwrite("splice.gff", bed);
  }

  for(i=0; i < map.noofsplits; i++) {
    cidx = map.cidx[i];
    chr = bl_fastaGetDescription(fasta, cidx);
    pos = map.pos[i];
 
    len = strlen(chr);
    chr2 = ALLOCMEMORY(space, NULL, char, len+1);
    memmove(chr2, chr, len);
    chr2[len] = 0;
    ptr = strtok(chr2, " ");

    //printsplits(space, ptr, pos, &map.cs[i], &map);
    FREEMEMORY(space, chr2);
  }

  bl_matchfileDestructSplitMap(space, &map);
  bl_matchfileDestructSpliceMap (space, sm);  
  FREEMEMORY(space, sm);

  return ;
}

/*----------------------------------- main -----------------------------------
 *    
 * @brief the main
 * @author Steve Hoffmann 
 *   
 */
 

int
main(int argc, char **argv) {
  
  void *space = NULL;
  
  manopt_optionset optset;
  manopt_arg *unflagged; 
  manopt_arg *dbfilenames;
  manopt_arg *queries;
  manopt_arg *indices;
  manopt_arg *grouplist;
  manopt_arg *bedfilenames;

  matchfile_t **files = NULL;
  annotationtrack_t *track = NULL;
  fasta_t *fasta = NULL;
  matchfileindex_t *idx2 = NULL;
  Uint prefixlen=0; 
  Uint splicebasenamelen;
  Uint minsplitno = 4;
  unsigned char gzip = 0, browse=0, test=0, call=0, saveindex=0, nostats=0, dumpstats=0, starts=0, wiggle=0, splice=0;
  char version[]="0.1";
  char *splicebasename = NULL;
  int *groups = NULL;
  int i;

  FILE *beddev = NULL;
  char *filename;
  char *transfilename;
   
  initIUPAC(1,1); 
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de"); 
  manopt(&optset, LISTOPT, 1, 'd', "database", 
      "list of path/filename(s) of database sequence(s)", "<file> [<file> ...]", 
      NULL, NULL);
  manopt(&optset, LISTOPT, 1, 'q', "query", 
      "path/filename of alignment file", "<file> [<file> ...]", NULL, NULL); 
  manopt(&optset, LISTOPT, 0, 'i', "index", 
      "path/filename of db index", "[<file> ... ]", NULL, NULL);
  manopt(&optset, LISTOPT, 0, 'a', "annotation", 
      "path/filename of bed annotation", "[<bedfile> ... ]", NULL, NULL);
  manopt(&optset, LISTOPT, 0, 'x', "generate", 
      "generate db index and store to disk", "<file>", NULL, NULL); 
  manopt(&optset, LISTOPT, 0, 'g', "group", 
      "group number for all nput files (1-offset and at most #files groups)", "[<int> ...]", NULL, NULL); 
  manopt(&optset, REQSTRINGOPT, 0, 's', "splice", 
      "dump splice sites to <basename>", NULL, NULL, &splicebasename);
  manopt(&optset, FLAG, 0, 'S', "starts", 
      "dump start sites", NULL, NULL, &starts);
  manopt(&optset, FLAG, 0, 'b', "browse", 
      "start browser", NULL, NULL, &browse);
  manopt(&optset, FLAG, 0, 'c', "call",
      "variant caller", NULL, NULL, &call);
  manopt(&optset, FLAG, 0, 'w', "expression",
      "generate a expression graph file from the match files", NULL, NULL, &wiggle);
  manopt(&optset, FLAG, 0, 'Z', "dumpstats",
      "dump data stats", NULL, NULL, &dumpstats);
  manopt(&optset, FLAG, 0, 'H', "nostats", 
      "no stats calculation", NULL, NULL, &nostats);
  manopt(&optset, FLAG, 0, 'T', "test",
      "test the program", NULL, NULL, &test);
  manopt(&optset, REQUINTOPT, 0, 'm', "minsplitno", 
      "minimum number of splits required to call a splice site", "<n>",
      NULL, &minsplitno);


  unflagged = manopt_getopts(&optset, argc, argv);
  saveindex = manopt_isset(&optset, 'x', NULL);
  
  if(!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
    manopt_help(&optset, "please give index filename using -i XOR -x option\n");
  } else if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unknown argument(s)\n");
  }

  MSG("reading database sequences.\n"); 
  NFO("minsplitno set to %d\n", minsplitno);

  dbfilenames = manopt_getarg(&optset, 'd', "database");
  fasta = bl_fastxGetSet(space, dbfilenames->values, 
      dbfilenames->noofvalues, 1, 0, 0, 1);

  NFO("%d database sequences found.\n", fasta->noofseqs);
  MSG("reading query files.\n");

  queries = manopt_getarg(&optset, 'q', "query");
  if(queries->noofvalues > 30) {
    manopt_help(&optset, "currently no more than 30 query files allowed\n");
  }

  grouplist = manopt_getarg(&optset, 'g', "group");
  if(grouplist) {
    if(grouplist->noofvalues != queries->noofvalues) 
      manopt_help(&optset, "please provide a group name for each input file");

    groups = ALLOCMEMORY(space, NULL, Uint, grouplist->noofvalues);
    
    for(i=0; i < grouplist->noofvalues; i++) {
      groups[i] = atoi(grouplist->values[i]);
      if(groups[i] == 0)
        manopt_help(&optset, "please provide group numbers (int) > 0");
      if(groups[i] > queries->noofvalues)
        manopt_help(&optset, "please provide groupnumbers <= number of input files");
      NFO("found group number %d\n", groups[i]);
    }
  }

  if(saveindex) {
    indices = manopt_getarg(&optset, 'x', "generate");
  } else {
    indices = manopt_getarg(&optset, 'i', "index");
  }

  if(indices->noofvalues != queries->noofvalues) {
    manopt_help(&optset, "please provide an index file name for each query file\n");
  }

  ntcode  = getNTcodekey(space);
  files = ALLOCMEMORY(space, NULL, matchfile_t*, queries->noofvalues);

  for(i=0; i < queries->noofvalues; i++) {

    files[i] = ALLOCMEMORY(space, NULL, matchfile_t, 1);  
    files[i]->fmt = 0;
    files[i]->index = NULL;
    files[i]->filename = queries->values[i];

    prefixlen = bl_fileprefixlen(files[i]->filename);

    if(strncmp(&files[i]->filename[prefixlen], ".gz", 3) == 0 || 
        strncmp(&files[i]->filename[prefixlen], ".gzip", 3) == 0) {
      gzip = 1;
    }

    files[i]->gzip = gzip;

    if(saveindex) {

      bl_matchfileIndex(space, files[i], fasta);
      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
      if(!nostats) { 
        if(!test) { 
          getstats(space, files[i], fasta); 
        } else { 
          teststats(space, files[i],fasta);
        }
      }

      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
      idx2 = bl_matchfileReadIndex(space, indices->values[i]);

      fprintf(stderr, "compare index (%p:%p):%d\n", 
          (void*)files[i]->index, (void *)idx2, 
          bl_compareIndices(files[i]->index, idx2));
      bl_matchfileDestructIndex(space, idx2); 
      FREEMEMORY(space, idx2);

    } else if(indices->values[i]) {
    
      MSG("reading index file\n");
      files[i]->index = bl_matchfileReadIndex(space, indices->values[i]);
    }
   
    /*
     *  typically stats should be present if stats
     *  failed however we are not allowed to set
     *  this one
     */
    if(files[i]->index->stats)
    files[i]->index->stats->mincover = 6;
  }

  if(manopt_isset(&optset, 'a', "annotation")) { 
    bedfilenames = manopt_getarg(&optset, 'a', "annotation");
    for(i=0; i < bedfilenames->noofvalues; i++) { 
  
      prefixlen = bl_fileprefixlen(bedfilenames->values[i]);
      if( strncmp(&bedfilenames->values[i][prefixlen], ".bed", 3) == 0 ||
        strncmp(&bedfilenames->values[i][prefixlen], ".BED", 3) == 0) {
        track = bl_BEDread(space, bedfilenames->values[i]);

        beddev = fopen("sorted.bed", "w");
        if(beddev  == NULL) {
          fprintf(stderr, "could not open file %s. Exit forced.", "sorted.bed");
          exit(-1);
        }

        bl_BEDwrite(track, beddev);
      
      } else if( strncmp(&bedfilenames->values[i][prefixlen], ".gff", 3) == 0 ||
        strncmp(&bedfilenames->values[i][prefixlen], ".GFF", 3) == 0) {
        track = bl_GFFread(space, bedfilenames->values[i]);
        bl_BEDwrite(track, beddev);
      
      } else {
        manopt_help(&optset, "please provide files with .GFF or .BED extension\n");
      }
    }
  }

  /*
   * TODO: stats for all files, what if stats failed!
   */

  if(dumpstats) bl_matchfileDumpSampleStats(files[0]->index->stats);
  if(call) eval(space, files, groups, queries->noofvalues, fasta);

  if(manopt_isset(&optset, 's', "splice") && 
      manopt_isset(&optset, 'q', "query")) {

      if(!splicebasename) { 
        splicebasename = bl_basename(files[0]->filename);
      }

      splicebasenamelen  = strlen(splicebasename);

      filename = ALLOCMEMORY(space, NULL, char, splicebasenamelen+7+4+1);
      sprintf(filename, "%s.splice.bed", splicebasename);

      transfilename = ALLOCMEMORY(space, NULL, char, splicebasenamelen+7+4+1);
      sprintf(transfilename, "%s.trans.bed", splicebasename);

      splice = 1;
      evalsplice(space, files, groups, queries->noofvalues, fasta, track, filename, transfilename, minsplitno);
  }

  
  if(browse) view(space, files, queries->noofvalues, fasta, track);
  if(starts) evalstarts(space, files, groups, queries->noofvalues, fasta, track);
  if(wiggle) writeexpression(space, files, groups, queries->noofvalues, fasta, track);

  bl_fastaDestruct(space, fasta);
  FREEMEMORY(space, fasta);
  
  if(files) {
    for(i=0; i < queries->noofvalues; i++) { 
      if (files[i]->index) {
        bl_matchfileDestructIndex(space, files[i]->index);
        FREEMEMORY(space, files[i]->index);
      }
      FREEMEMORY(space, files[i]);
    }
    FREEMEMORY(space, files);
  }

  if(groups) {
    FREEMEMORY(space, groups);
  }


  if(track) { 
    bl_annotationtrackDestruct(space, track);
    FREEMEMORY(space, track);
  }


  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(space, unflagged);

  return 0;
}

