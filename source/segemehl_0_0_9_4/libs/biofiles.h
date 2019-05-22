#ifndef _BIOFILES_
#define _BIOFILES_

/*
 *
 *	biofiles.h
 *  declarations
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/10/2007 02:32:29 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: biofiles.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/biofiles.h $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"

typedef struct {

  CharSequence** seqs;
  Uint noofseqs;
  Uint minlen;
  Uint maxlen;

} fasta_t;

typedef struct {
    char *seqname;
    Uint seqnamelen;
    char *source;
    Uint sourcelen;
    char *feat;
    Uint featlen;
    Uint start;
    Uint end;
    double score;
    unsigned char strand;
    unsigned char frame;
    char* attrib;
    Uint attriblen;

} gff_t;

void destructfasta(void *space, fasta_t* f);
fasta_t** chopfasta(void *space, fasta_t* f, Uint pieces);
fasta_t* initfasta(void *);
void addfasta(void *space, fasta_t*, char *desc, Uint, char* sequence, Uint);
fasta_t* readfasta(void *space, char* filename, unsigned char upper, unsigned char lower);
fasta_t* readsolexa(void *space, char* filename);
fasta_t*getfastaset(void *space, char **filenames, unsigned int nooffiles,
    unsigned char upper, unsigned char lower);
 
#endif
