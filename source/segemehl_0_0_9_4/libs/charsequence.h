 #ifndef INTSEQUENCE_H
 #define INTSEQUENCE_H

/*
 * charsequence.h
 * declaration of char sequence
 * and functions working on it
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 87 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-20 11:24:26 +0100 (Thu, 20 Nov 2008) $
 *
 *  Id: $Id: charsequence.h 87 2008-11-20 10:24:26Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/charsequence.h $
 */

 #include "basic-types.h"
 #include <stdio.h>
 #include <stdlib.h>

 typedef struct {
	Uint descrlen;			
    Uint namelen;
	Uint urllen;
    /*new: num of nfo*/
    Uint noofinfo;
	
	char *description; 		/*a description*/
    char *alphabetname;		/*the name of the corresponding alphabet*/
	char *url;				/*the name of the sequences url*/
	
	char *sequence;			/*the sequence itself*/
	/*old: Uint *info*/
    char **info;		    /*additional information*/
	Uint length;

    Uint *map;
    Uint mapsize;
	
 } CharSequence;
 
 void destructSequence(void *, CharSequence *);
 CharSequence* initSequence(void *);
 void resetSequence(CharSequence *);
 char* printSequence(void *, CharSequence *, Uint);
 void dumpSequence(CharSequence *s);
 void saveSequence (CharSequence *s, char *filename); 
 CharSequence* loadSequence (void *space, char *filename);
 char * printAlignment (void *, int *, Uint, CharSequence *, CharSequence *, 
	 Uint);
 CharSequence** createSequenceHash(void *, Uint);

 static inline char* charDNAcomplement(void *space, char *s, Uint len) {
    Uint i,k=0;
    char* buffer;
    
    buffer = ALLOCMEMORY(space, NULL, char, len+1);
    for(i=len; i > 0; i--) {
        switch(s[i-1]) {
          case 'a':
            buffer[k] = 't';
            break;
          case 't':
            buffer[k] = 'a';
            break;
          case 'c':
            buffer[k] = 'g';
            break;
          case 'g':
            buffer[k] = 'c';
            break;
          case 'n':
            buffer[k] = 'n';
            break;
          case 'A':
            buffer[k] = 'T';
            break;
          case 'T':
            buffer[k] = 'A';
            break;
          case 'C':
            buffer[k] = 'G';
            break;
          case 'G':
            buffer[k] = 'C';
            break;
          case 'N':
            buffer[k] = 'N';
            break;
          default:
            buffer[k] = s[i-1]; 
            break;
        }
     k++;
    }
    buffer[k] = '\0';
    return buffer;
 }
 
 #endif
