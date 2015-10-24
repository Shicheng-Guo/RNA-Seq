
/*
 *  biofiles.c
 *  helper functions to handle file types
 *  used in bioinformatics
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 01:56:15 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 76 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-11 16:34:21 +0100 (Tue, 11 Nov 2008) $
 *
 *  Id: $Id: biofiles.c 76 2008-11-11 15:34:21Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/biofiles.c $
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileio.h"
#include "charsequence.h"
#include "assert.h"

/*-------------------------------- readfasta ---------------------------------
 *    
 * @brief reads a fasta file 
 * @returns a stringset containing seqs in fastafile
 * @author Steve Hoffmann 
 *   
 */


fasta_t* 
initfasta(void *space) {
  fasta_t *f;

  f = ALLOCMEMORY(space, NULL, fasta_t, 1);
  f->seqs = NULL;
  f->noofseqs = 0;
  f->minlen = 0;
  f->maxlen = 0;

  return f;
}

void
destructfasta(void *space, fasta_t* f) {
  Uint i;
    
  for(i=0; i < f->noofseqs; i++) {
    destructSequence(space, f->seqs[i]);
  }

  FREEMEMORY(space, f->seqs);
}

void
addfasta(void *space,
    fasta_t *f,
    char    *desc,
    Uint    descrlen,
    char    *sequence,
    Uint    seqlen) {


  f->seqs = ALLOCMEMORY(space, (f->seqs), CharSequence*, (f->noofseqs)+1);
  f->seqs[f->noofseqs] = initSequence(space);
  f->seqs[f->noofseqs]->description = desc;
  f->seqs[f->noofseqs]->descrlen = descrlen;
  f->seqs[f->noofseqs]->sequence = sequence;
  f->seqs[f->noofseqs]->length = seqlen;
  f->minlen = (seqlen < f->minlen) ? seqlen : f->minlen;
  f->maxlen = (seqlen > f->maxlen) ? seqlen : f->maxlen;

  f->noofseqs++;
  return;
}

fasta_t**
chopfasta(void *space, fasta_t* f, Uint pieces) {
  Uint size, r, i, j, offset=0;
  fasta_t **chops; 

  size = f->noofseqs/pieces;
  r = f->noofseqs-(pieces*size);

  assert((pieces*size)+r == f->noofseqs);

  chops = ALLOCMEMORY(space, NULL, fasta_t*, pieces);
  for(i=0; i < pieces; i++) {

    chops[i] = initfasta(NULL);
    if (i < pieces-1) 
      chops[i]->noofseqs=size;
    else
      chops[i]->noofseqs = size+r;

    chops[i]->seqs = ALLOCMEMORY(space, NULL, CharSequence*, chops[i]->noofseqs);

    for(j=0; j < chops[i]->noofseqs; j++) {  
      chops[i]->seqs[j] = f->seqs[j+offset];
      chops[i]->minlen = (f->seqs[j+offset]->length < chops[i]->minlen) ? 
        f->seqs[j+offset]->length : chops[i]->minlen;
      chops[i]->maxlen = (f->seqs[j+offset]->length > chops[i]->maxlen) ? 
        f->seqs[j+offset]->length : chops[i]->maxlen;
    }

    offset += chops[i]->noofseqs;
  }
  return chops;
}


fasta_t* 
readfasta(void *space, char* filename, unsigned char upper, unsigned char lower) {

  char ch;
  char *buffer;
  char *descrbuffer = NULL;
  Uint  descrlength = 0;
  FILE *fp;
  unsigned char desc=1;
  unsigned char first=0;
  fasta_t *fasta;
  Uint buffersize = MAXBUFFERSIZE;
  Uint len=0;  

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  fasta = initfasta(space);

  fp = fopen(filename, "r");

  if (fp == NULL) {
    fprintf(stderr,"Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }
  
  while((ch=getc(fp)) != EOF) {	

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch=='>' && !first) {
        desc=1;
        first = 1;
    }

    if(ch=='>' && len > 0) {
      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len]='\0';
      desc = 1;
      addfasta(space, fasta, descrbuffer, descrlength, buffer, len);

      len = 0;
      descrlength = 0;
      descrbuffer = NULL;
      buffersize = MAXBUFFERSIZE;
      buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    }

    if(!desc && ch =='\n') {
      /*do nothing.*/
    } else {
      if(desc && ch == '\n') { 
        buffer = ALLOCMEMORY(space, buffer, char, len+1);  
        buffer[len]='\0'; 

        descrbuffer = buffer;
        descrlength = len;

        len = 0;
        buffersize = MAXBUFFERSIZE;
        buffer = ALLOCMEMORY(space, NULL, char, buffersize);
        desc = 0;
      } else {
        len++;

        if (upper && !desc) {   
          buffer[len-1]=(char)toupper((int)ch);
        } else if (lower && !desc) { 
          buffer[len-1]=(char)tolower((int)ch);
        } else { 
          buffer[len-1]=(char) ch;
        }

      }
    }
  }

  buffer = ALLOCMEMORY(space, buffer, char, len+1);  
  buffer[len]='\0'; 
  addfasta(space, fasta, descrbuffer, descrlength, buffer, len);

  fclose(fp);
  return fasta;
}


fasta_t*
getfastaset(void *space, char **filenames, unsigned int nooffiles,
    unsigned char upper, unsigned char lower) {
  
  unsigned int i;
  fasta_t *file, *set;

  set = initfasta(space);
  for(i=0; i < nooffiles; i++) {
    file = readfasta(space, filenames[i], upper, lower);
    set->seqs = ALLOCMEMORY(space, set->seqs, CharSequence*, 
        set->noofseqs+file->noofseqs);
    memmove(&set->seqs[set->noofseqs], file->seqs, 
        file->noofseqs*sizeof(CharSequence*));
    set->noofseqs += file->noofseqs;
    FREEMEMORY(space, file->seqs);
    FREEMEMORY(space, file);
  }
  return set;
}


gff_t*
initGff(void *space) {
  gff_t *g;

  g = ALLOCMEMORY(space, NULL, gff_t, 1);
  g->seqname = NULL;
  g->seqnamelen = 0;
  g->source = NULL;
  g->sourcelen = 0;
  g->feat = NULL;
  g->featlen = 0;
  g->start = 0;
  g->end = 0;
  g->score = .0;
  g->strand = '0';
  g->frame = '0';
  g->attrib = NULL;
  g->attriblen = 0;

  return g;
}


void
writeGff(char *filename, gff_t *set, Uint len) {
  Uint i;
  gff_t *g;
  FILE *fp;


  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    g = &set[i];  
    fprintf(fp,"%s\t%s\t%s\t", g->seqname, g->source, g->feat);
    fprintf(fp,"%d\t%d\t%c\t", g->start, g->end, g->strand);
    fprintf(fp,"%c\t%s\n", g->frame, g->attrib);
  }

  fclose(fp);
  return;
}



