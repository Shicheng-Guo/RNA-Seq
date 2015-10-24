#ifndef EVALMATCHFILE_H
#define EVALMATCHFILE_H

/*
 *
 *	evalmatchfiles.h
 *  evalutate matchfiles
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/14/2010 12:07:57 AM CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "matchfiles.h"
#include "biofiles.h"


#define MAX_STARTSITES 20000
#define EPSILON_NV 0.0000001
#define SDFRACTION 3
#define MINSUBPROB -4.0


typedef struct {
  Uint noofgroups;
  char *type;
  double *s_cons;
  double *s_ref;
  double *s_consx;
  double *s_refx;
  char **chars;
  Uint *len;
} matchfileTestGroups_t;

extern char * ntcode;
extern char * ntdecode;

  void
bl_matchfileEvalCrossSections (void *space,  matchfile_t **file, int *gropus, Uint nooffiles, fasta_t *set, 
    Uint (*f)(void *, Uint fidx, Uint cidx, Uint pos, matchfileCross_t*, char, matchfileSampleStats_t *, unsigned char, void *), void *nfo);

void bl_matchfileGetConsensus(matchfileFrame_t *frame);
Uint* bl_matchfileGetNTCounts(matchfileCross_t *cs);
double* bl_matchfileGetNTError(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTRedundancy(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTEdist(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTReadPos(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTReadPosVar(matchfileCross_t *);
matchfileFrameStats_t* bl_matchfileFrameStats(void *space, matchfileFrame_t *frame);
void bl_matchfileDestructFrameStats(void *space, matchfileFrameStats_t *stats);
void bl_matchfileRSSGNUPLOT(void *space, matchfileFrame_t *frame, matchfileFrameStats_t *stats);
void bl_matchfileCOVGNUPLOT(void *space, matchfileFrame_t *frame);
extern FILE *popen( const char *command, const char *modes);
extern int pclose(FILE *stream);
Uint bl_matchfileSampleCmp (Uint elemA, Uint elemB, void *toSort, void *info);
void bl_matchfileGetErrorDensity(void *space, matchfileFrame_t *frame, Uint pos, matchfileFrameStats_t *, void *nfo);
Uint bl_matchfileTest(void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, char ref, matchfileSampleStats_t *, unsigned char, void *nfo);
int bl_matchfileSampleCrossSections(void *space, matchfile_t *file, fasta_t *set, Uint n, 
    void (*f)(void *, matchfileFrame_t*, Uint, matchfileFrameStats_t *, void *), void *info);
void bl_matchfileGetConditionals (void *space, matchfileFrame_t *frame,
    Uint pos, matchfileFrameStats_t *stats, void *nfo);
matchfileSampleStats_t*
bl_matchfileInitSampleStats (void *space, Uint maxsample, Uint maxcover, Uint mincover, 
    Uint areasize, double maxareae);
void bl_matchfileDestructSampleStats (void *space, matchfileSampleStats_t *stats);
void bl_matchfileDumpSampleStats (matchfileSampleStats_t *stats);
#endif
