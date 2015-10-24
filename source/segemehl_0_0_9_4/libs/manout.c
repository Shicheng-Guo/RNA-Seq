/*
 * manout.c
 * attempt for flexible output of genome mapping w/ SEGEMEHL
 *
 * @author Christian Otto
 * @email christan@bioinf.uni-leipzig.de
 * @date Wed Sep 24 10:56:23 CEST 2008
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "basic-types.h"
#include "bitArray.h"
#include "memory.h"
#include "mathematics.h"
#include "sort.h"
#include "info.h"
#include "biofiles.h"
#include "fileio.h"
#include "vtprogressbar.h"
#include "debug.h"
#include "charsequence.h"
#include "manout.h"
#include <assert.h>
#include "alignment.h"
#include "manoutformats.h"

/*-------------------------------- initStruct ---------------------------------
 *    
 * @brief inits
 * @author Steve Hoffmann 
 *   
 */
 
inline void
initMatch (gmatch_t* m)
{
    m->i = 0; 
    m->p = 0;
    m->q = 0;
    m->checklen = 0;
    m->scr = 0;
    m->mat = 0;
    m->mis = 0;
    m->del = 0;
    m->ins = 0;
    m->subject = 0;
    m->edist = 0;
    m->al = NULL;
    return ;
}

void 
initRead(gread_t *r, CharSequence *query) 
{
    r->query = query;
    r->noofmatches_plus = 0;
    r->noofmatches_minus = 0;
    r->plus = NULL;
    r->minus = NULL;

    return;
}

void
initGmap(Gmap *map, MultiCharSeq *mseq, Uint offset) 
{ 
  map->mseq = mseq;
  map->mapoffset = offset;
  map->noofreads = 0;
  map->reads = 0;

  return;
}

inline void
setMatches(gread_t *r, gmatch_t *m, Uint noofmatches,
    unsigned char strand){

  if (strand == MINUSSTRAND)
  {
    r->noofmatches_minus = noofmatches;
    r->minus = m;
  } else {
    r->noofmatches_plus = noofmatches;
    r->plus = m;
  }

  return;
}



inline void
setReads(Gmap *map, gread_t *r, Uint noofreads){
 map->reads = r;
 map->noofreads = noofreads;
}


/*------------------------------- matchHeader --------------------------------
 *    
 * @brief write match header
 * @author Steve Hoffmann 
 *   
 */
 
void
matchHeader (FILE *dev, Uint rep_type)
{
  fprintf(dev, HEAD[rep_type]);
  return ;
}


/*----------------------------- genericOutput --------------------------------
 *    
 * @brief write generic output with format string as parameter
 * @author Christian Otto
 *   
 */

void
genericOutput (FILE *dev, const char *format, ...){
  va_list args;
  va_start(args, format);
  vfprintf(dev, format, args);
  va_end(args);
}


/*------------------------------- reportMatch --------------------------------
 *    
 * @brief reports a match to device with different output formats
 * @author Steve Hoffmann 
 *   
 */

inline void
reportMatch (FILE* dev, Gmap *map, Uint rep_type, unsigned char showalignment, pthread_mutex_t* mtx, Uint pmaxedist, Uint mmaxedist){

  int i, j, seqlen, qrylen, size;
  Uint off, seqstart=0;
  gread_t *reads,
          *read;
  char *seq,
       *qry,
    *evalue,
    *meopstr;

  const char *fmt;

  if (mtx != NULL)  {
    pthread_mutex_lock(mtx);
  }
  
  fmt = FORMAT[rep_type];
  off = map->mapoffset;
  reads = map->reads;
  
  for(i = 0; i < map->noofreads; i++) {
    read = &reads[i];

    for(j = 0; j < read->noofmatches_plus; j++) {
      if(read->plus[j].subject > 0) {
        seqstart = map->mseq->markpos[read->plus[j].subject-1]+1; 
      } else {
        seqstart = 0;
      }
      if(read->plus[j].edist > pmaxedist) {
        if (read->plus[j].al) {
            wrapAlignment(read->plus[j].al);
            FREEMEMORY(space, read->plus[j].al);
        }
        continue;
      }
      if (off+read->plus[j].checklen < read->plus[j].j) {
	DBGL(1, "failed: checklen: %s\n", read->query->description);  
      }
      qrylen = read->query->length+1;
      qry = (char *) malloc((qrylen + 1) * sizeof(char));
      memset(qry, (char) 0, qrylen + 1);
      memmove(qry, read->query->sequence, qrylen);
      
      seqlen = off+read->plus[j].q - read->plus[j].p;
      seq = (char *) malloc((seqlen + 1) * sizeof(char));
      memset(seq, (char) 0, seqlen + 1);
      memmove(seq, map->mseq->sequences + read->plus[j].p, seqlen);
      size = snprintf(NULL, 0, "%.4f", read->plus[j].evalue);
      evalue = (char *) malloc((size + 1) * sizeof(char));
      snprintf(evalue, size + 1, "%.4f", read->plus[j].evalue);      
      //use alignment instead
      if (seqstart > off+read->plus[j].p) continue;
      meopstr = multieopstring(read->plus[j].al);

#ifdef DIRECTOUT

fprintf(dev, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\t%s\n", 
                         read->query->description, 
                         read->plus[j].edist, 
                         read->plus[j].scr,
                         evalue,
                         off+read->plus[j].i,
                         off+read->plus[j].j,
                         read->plus[j].mat,
                         read->plus[j].mis,
                         read->plus[j].ins,
                         read->plus[j].del,
                         '+',
                         off+read->plus[j].p-seqstart, 
                         off+read->plus[j].q-seqstart,
                         ((CharSequence*)map->mseq->ref[read->plus[j].subject].ref)->description, 
                         meopstr);
#else
      
      genericOutput(dev, fmt, read->query->description, read->query->length, 
		    read->plus[j].scr, evalue,
		    off+read->plus[j].i, off+read->plus[j].j,
		    off+read->plus[j].p-seqstart, off+read->plus[j].q-seqstart,
		    read->plus[j].mat, read->plus[j].mis,
		    read->plus[j].ins, read->plus[j].del, '+',
		    read->plus[j].edist,
            seq, 
            ((CharSequence*)map->mseq->ref[read->plus[j].subject].ref)->description, 
            qry, meopstr);
#endif
      free(meopstr);
      free(seq);
      free(qry);
      free(evalue);
      if(read->plus[j].al) {
        if (showalignment)
        showAlign(read->plus[j].al, dev);
        wrapAlignment(read->plus[j].al);
        FREEMEMORY(space, read->plus[j].al);
      }
    }
    for(j = 0; j < read->noofmatches_minus; j++) {
       
      if(read->minus[j].subject > 0) {
        seqstart = map->mseq->markpos[read->minus[j].subject-1]+1; 
      } else {
        seqstart = 0;
      }
      
      if(read->minus[j].edist > mmaxedist) {
        if (read->minus[j].al) {
           wrapAlignment(read->minus[j].al);
           FREEMEMORY(space, read->minus[j].al);
        }
        continue;
      }

      if (off+read->minus[j].checklen < read->minus[j].j) {
	DBGL(1, "failed: checklen: %s\n", read->query->description);  
      }
      qrylen = read->query->length+1;
      qry = (char *) malloc((qrylen + 1) * sizeof(char));
      memset(qry, (char) 0, qrylen + 1);
      memmove(qry, read->query->sequence, qrylen);
      
      seqlen = off+read->minus[j].q - read->minus[j].p;
      seq = (char *) malloc((seqlen + 1) * sizeof(char));
      memset(seq, (char) 0, seqlen + 1);
      memmove(seq, map->mseq->sequences + read->minus[j].p, seqlen);
      size = snprintf(NULL, 0, "%.4f", read->minus[j].evalue);
      evalue = (char *) malloc((size + 1) * sizeof(char));
      snprintf(evalue, size + 1, "%.4f", read->minus[j].evalue);
      //use alignment instead
      if (seqstart > off+read->minus[j].p) continue; 
      meopstr = multieopstring(read->minus[j].al);

#ifdef DIRECTOUT

fprintf(dev, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\t%s\n", 
                         read->query->description, 
                         read->minus[j].edist, 
                         read->minus[j].scr,
                         evalue,
                         off+read->minus[j].i,
                         off+read->minus[j].j,
                         read->minus[j].mat,
                         read->minus[j].mis,
                         read->minus[j].ins,
                         read->minus[j].del,
                         '-',
                         off+read->minus[j].p-seqstart, 
                         off+read->minus[j].q-seqstart,
                         ((CharSequence*)map->mseq->ref[read->minus[j].subject].ref)->description, 
                         meopstr);
#else
      
      genericOutput(dev, fmt, read->query->description, read->query->length,
		    read->minus[j].scr, evalue,
		    off+read->minus[j].i, off+read->minus[j].j,
		    off+read->minus[j].p-seqstart, off+read->minus[j].q-seqstart, 
		    read->minus[j].mat, read->minus[j].mis,
		    read->minus[j].ins, read->minus[j].del, '-',
		    read->minus[j].edist,
            seq, 
            ((CharSequence*)map->mseq->ref[read->minus[j].subject].ref)->description, 
            qry, meopstr);
#endif
      free(meopstr);
      free(seq);
      free(qry);
      free(evalue);
      if(read->minus[j].al) {
        if (showalignment)
        showAlign(read->minus[j].al, dev);
        wrapAlignment(read->minus[j].al);
        FREEMEMORY(space, read->minus[j].al);
      }
    }
  }

  if (mtx != NULL){
    pthread_mutex_unlock(mtx);
  }

  return ;
}

