#ifndef MERGE_H
#define MERGE_H

/*
 * merge.h
 * functions to merge matches
 *
 *  SVN
 *  Revision of last commit: $Rev: 310 $
 *  Author: $Author: steve $
 *  Date: $Date: 2011-07-15 17:03:40 +0200 (Fr, 15. Jul 2011) $
 *
 *  Id: $Id: merge.h 310 2011-07-15 15:03:40Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/merge.h $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fileBins.h"
#include "biofiles.h"

#define SAM 0

typedef struct {
  /* match (always first one read,
   * not necessarily 1st in pair) */
  char *match;
  /* match length */
  Uint matchlen;
  /* mate match */
  char *matematch;
  /* mate match length */
  Uint matematchlen;
  /* key for merging */
  char *key;
  /* key length */
  Uint keylen;
  /* flag in SAM format */
  Uint flag;
  /* mate flag in SAM format */
  Uint mateflag;
  /* edit distance of qry alignment */
  Uint edist;
  /* edit distance of mate alignment */
  Uint mateedist;
  /* number of matches */
  Uint noofmatches;
  /* number of mate matches */
  Uint noofmatematches;
  /* reference name */
  char *rname;
  /* reference start */
  Uint rstart;
  /* strand */
  char strand;
  /* mate reference name */
  char *matername;
  /* mate reference start */
  Uint materstart;
  /* mate strand */
  char matestrand;
} bl_mergefilematch_t;

typedef struct {
  /* file pointer */
  FILE *fp;
  /* EOF read */
  unsigned char eof;
  /* current entry read */
  bl_mergefilematch_t *entry;
  /* current entry complete? */
  unsigned char complete;
} bl_mergefile_t;

typedef struct {
  Uint nooffiles;
  bl_mergefile_t *f;
} bl_mergefiles_t;

void bl_mergefilesInit(void *space, bl_mergefiles_t *files, Uint nooffiles);
void bl_mergefilesDestruct(void *space, bl_mergefiles_t *files);
void bl_mergefileInit(void *space, bl_mergefile_t *file, FILE *fp);
void bl_mergefileDestruct(void *space, bl_mergefile_t *file);
void bl_mergefilematchInit(void *space, bl_mergefilematch_t *match);
int bl_mergefilematchCompare(bl_mergefilematch_t *i, bl_mergefilematch_t *j);
void bl_mergefilematchDestruct(void *space, bl_mergefilematch_t *match);
unsigned char bl_mergeParseLine(void *space, bl_mergefilematch_t *match, char *line, Uint len);
void bl_mergeReadNext(void *space, bl_mergefile_t *file);
void bl_mergeUpdateTag(void *space, char **line, Uint *len, Uint noofmatches);
void se_mergeBisulfiteBins (void *space, bl_fileBinDomains_t *bsdomains, fasta_t **reads,
			    FILE *dev, bl_fileBinDomains_t *chrdomains, unsigned char remove);
#endif /* MERGE_H */
