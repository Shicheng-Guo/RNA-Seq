#ifndef ALIGNMENT_H
#define ALIGNMENT_H

/*
 *
 *	alignment.h
 *  alignment representation
 *  
 *  idea: 
 *  Stephan Kurtz, Gordon Gremme. Foundations
 *  of sequence analysis, University Hamburg, 2005
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 02/03/2009 11:56:27 AM CET  
 *
 */

#include "basic-types.h"

typedef enum 
{
  Replacement, Deletion, Insertion
} Eoptype;

typedef struct
{
  Eoptype eop;
  Uint steps;
} Multieop;

typedef struct {
  char *u;
  char *v;
  Uint ulen;
  Uint vlen;

  /*start of aligment (use in approx string matching, local align)*/
  Uint uoff;
  Uint voff;
  Multieop *meops;
  Uint numofmeops;

} Alignment;


void showmultieoplist(Alignment *al);
void showDynmultieoplist(Alignment* al, int size);
void showAlign(Alignment* al, FILE *dev);
void initAlignment(Alignment *al, char *u, Uint ulen, Uint uoff, char *v, Uint vlen, Uint voff);
void insertEop(Alignment *al, Eoptype eop);
void revMeops(Alignment *al);
void wrapAlignment(Alignment *al);
Uint getEdist(Alignment *al);
void countEops(Alignment *al, Uint *mat, Uint *mis, Uint *ins, Uint *del);
char * multieopstring(Alignment *al);
#endif
