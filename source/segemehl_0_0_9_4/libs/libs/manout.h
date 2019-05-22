#ifndef MANOUT_H
#define MANOUT_H

/*
 * manout.h
 * attempt for flexible output of genome mapping w/ SEGEMEHL
 *
 * @author Christian Otto
 * @email christan@bioinf.uni-leipzig.de
 * @date Wed Sep 24 10:56:23 CEST 2008
 *
 */

#include <pthread.h>
#include "basic-types.h"
#include "bitArray.h"
#include "biofiles.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "alignment.h"

#define MINUSSTRAND 0
#define PLUSSTRAND 1



typedef struct gmatch_s{
  
  Uint i;
  Uint j;
  Uint p;
  Uint q;
  Uint checklen;
 
  int scr; 
  double evalue;
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  Alignment *al;
  Uint subject;

} gmatch_t;


typedef struct gread_s{

  CharSequence *query;

  Uint noofmatches_plus;
  gmatch_t* plus;
  Uint noofmatches_minus;
  gmatch_t* minus;

} gread_t;


typedef struct Gmap_s{

  MultiCharSeq *mseq;
  Uint mapoffset;
  Uint noofreads;
  gread_t *reads;

} Gmap;

extern void reportMap(FILE*, Gmap *map, Uint level);
extern void initMatch(gmatch_t *);
void initRead(gread_t *, CharSequence *);
void initGmap(Gmap *, MultiCharSeq *, Uint);
extern void setMatches(gread_t*, gmatch_t *, Uint, unsigned char);
extern void setReads(Gmap *, gread_t *, Uint);
extern void reportMatch (FILE* dev, Gmap *map, Uint level, unsigned char showalignment, pthread_mutex_t *mtx, Uint mmaxedist, Uint pmaxedist);
void matchHeader(FILE* dev, Uint level);
void genericOutput(FILE* dev, const char *fmt, ...);

#endif /* MANOUT_H */
