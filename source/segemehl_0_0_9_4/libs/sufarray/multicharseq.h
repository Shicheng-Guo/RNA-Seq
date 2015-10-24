#ifndef MULTI_SEQ_H
#define MULTI_SEQ_H

/*
 *	multiseq.h
 *  declarations for a datastructure containing
 *  multiple integer sequences
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/11/06 15:09:15 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 65 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-22 00:14:54 +0200 (Mon, 22 Sep 2008) $
 *
 *  Id: $Id: multicharseq.h 65 2008-09-21 22:14:54Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/multicharseq.h $
 */


#include "basic-types.h" 
#include "charsequence.h"

#define MSEQ_BSEARCH_THRESHOLD 10

typedef struct {
  void *ref; 	
} SeqReference;


typedef struct {
  Uint numofsequences;
  Uint totallength;
  Uint *markpos;		/*markpos[i] is the position of a*/
  /*separator character between S_i and S_i+1*/
  char *sequences; 	/*array of concatenated sequences*/
  SeqReference *ref;  /*ref[i] points to the original sequence*/
  /*that starts at position markpos[i]*/
  char *map;
  Uint mapsize;

  char delim;

} MultiCharSeq;


void dumpMultiCharSeq (MultiCharSeq *);
MultiCharSeq* concatCharSequences(void *, CharSequence **, Uint, char, char);
void destructMultiCharSeq(void*, MultiCharSeq *);
Uint getMultiCharSeqIndex(MultiCharSeq *, char *);
Uint getMultiCharSeqRelPos(MultiCharSeq *, char *);
CharSequence* getCharSequence(MultiCharSeq *, Uint idx);

#endif
