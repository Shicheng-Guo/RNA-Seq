#ifndef SEGEMEHL_H
#define SEGEMEHL_H

/*
 *
 *	segemehl.h
 *  declarations for threaded segemehl
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 01/02/2008 10:12:46 AM CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 101 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-08 02:29:27 +0100 (Mon, 08 Dec 2008) $
 *
 *  Id: $Id: segemehl.h 101 2008-12-08 01:29:27Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.h $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include "biofiles.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include <pthread.h>

#define VERSION "  Thu Jan 14 13:11:37 CET 2010"

typedef struct segemehl_s {
    void *space;
    char *outfile;
    char *queryfilename;
    char *idxfilename;
    char *dbfilename;
    char *uninformativename;
    FILE *dev;
    FILE *uninformativedev;
    Suffixarray *sarray;
    fasta_t *reads;
    fasta_t *fasta;
    Uint totallength;
    Uint minsize;
    Uint *mapmatches;
    Uint *counter;
    Uint M;
    Uint s_ext;
    Uint p_mis;
    Uint Xoff;
    Uint k;
    Uint k_p;
    Uint rep_type;
    Uint check;
    Uint kMis;
    Uint threadno;
    Uint  bestonly;
    unsigned char matchingstat;
    unsigned char align;
    unsigned char mute;
    double maxevalue;
    int accuracy;
    Uint bedist;
    Uint fusion;
} segemehlinfo_t;

typedef struct checkthread_s {
    Uint noofseqs;
    Uint *counter;
} checkthreadinfo_t;

inline static void
setdefault(segemehlinfo_t *info) {
  info->dev = stdout;
  info->idxfilename = NULL;
  info->dbfilename = NULL;
  info->queryfilename = NULL;
  info->sarray = NULL;
  info->fasta = NULL;
  info->reads = NULL;
  info->outfile=NULL;
  info->totallength=0;
  info->counter=0;
  info->minsize = 12;
  info->k_p = 1;
  info->threadno = 1;
  info->accuracy = 85;
  info->M = 100;
  info->s_ext = 2;
  info->p_mis = 4;
  info->Xoff = 8;
  info->rep_type = 4;
  info->kMis = 0;
  info->mute = 0;
  info->matchingstat = 0;
  info->bestonly = 0;
  info->check = 0;
  info->maxevalue = 5;
  info->space = NULL;
  info->uninformativename = NULL;
  info->uninformativedev = NULL;
  info->bedist=0;
  info->align=0;
  info->fusion=0;
}

#endif
