#ifndef FILEIO_H
#define FILEIO_H

/*
 * fileio.h
 * declarations for file io
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.h $
 */

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#include "stringutils.h"
int bl_UnixSortMerge(void *space, char **filenames, Uint nooffiles, const char *fieldstring, char *outfile);
char * bl_getTempFile(char *tmp, Uint tmplen);
int bl_UnixSort(void *space, char *filename, const char *fieldstring);
char* readfile(void *, char *, Uint*);
stringset_t **readcsv(void *, char *, char*, Uint *);
void writeY(char *, double  *, Uint, Uint, Uint);
void writeXYUint(char *filename, Uint *X, Uint *Y, Uint len);
int bl_fgets(void *space, FILE *fp, char **str);
char * bl_basename (const char *name);
char* bl_replacealphanum(char *s, Uint len);
char*  bl_replacenonalphanum(char *s, Uint len);
int bl_fileprefixlen(char *filename);
int bl_UnixSortMerge(void *space, char **filenames, Uint nooffiles, const char *fieldstring, char *outfile);
void bl_writeFileHeader(char *filename, char *header) ;
void bl_freplace(char *filename, char oldchar, char newchar, char stop);
void bl_freplacearr(char *filename, char* oldchars, char *newchars, Uint len, char stop);
void bl_freplacestr(char *filename, char* str, Uint len, char stop);
void writeXYZ(char *filename, double *X, double *Y, double *Z, Uint len);
int bl_rm(void *space, char *filename);
void writeYUint(char *filename, Uint *Y, Uint len, Uint xoff, Uint yoff);
void writeYUintNorm(char *filename, Uint *Y, Uint len, Uint yoff);

#endif
