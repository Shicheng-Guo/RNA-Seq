/*
 * fileio.c
 * functions to manipulate and read files
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "fileio.h"

#ifndef DIR_SEPARATOR
#define DIR_SEPARATOR '/'
#endif

#if defined (_WIN32) || defined (__MSDOS__) || defined (__DJGPP__) || \
  defined (__OS2__)
#define HAVE_DOS_BASED_FILE_SYSTEM
#ifndef DIR_SEPARATOR_2 
#define DIR_SEPARATOR_2 '\\'
#endif
#endif

/* Define IS_DIR_SEPARATOR.  */
#ifndef DIR_SEPARATOR_2
# define IS_DIR_SEPARATOR(ch) ((ch) == DIR_SEPARATOR)
#else /* DIR_SEPARATOR_2 */
# define IS_DIR_SEPARATOR(ch) \
  (((ch) == DIR_SEPARATOR) || ((ch) == DIR_SEPARATOR_2))
#endif /* DIR_SEPARATOR_2 */



char* dirname(const char *filename) {
    char *s;

    s=strrchr(filename, (int)'/');
    if(s && *s)
      *s = '\0';

    return s;
}



  char *
basename (const char *name)
{
  const char *base;

#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
  /* Skip over the disk name in MSDOS pathnames. */
  if (ISALPHA (name[0]) && name[1] == ':') 
    name += 2;
#endif

  for (base = name; *name; name++)
  {
    if (IS_DIR_SEPARATOR (*name))
    {
      base = name + 1;
    }
  }
  return (char *) base;
}

char* readfile(void* space, char* filename, Uint* strlen) {

  char ch;
  char *buffer;
  FILE *fp;
  Uint buffersize = MAXBUFFERSIZE;
  Uint len=0;

  fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr, "Opening of file %s failed. Exit forced.\n", filename);
    exit(EXIT_FAILURE);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  while((ch=getc(fp)) != EOF) {
    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }
    len++;
    buffer[len-1]=(char)ch;	
  }
  buffer[len]='\0';
  fclose(fp);

  *strlen = len;
  return buffer;

}

stringset_t **
readcsv(void *space, 
    char* filename, 
    char *delim, 
    Uint *linecount) {

  Uint i, contentlen;
  char *content;
  stringset_t *lines, **csv;

  content = readfile(space, filename, &contentlen);
#ifndef _CRLF_
  lines = tokensToStringset(space, "\n", content, contentlen);
#else
  lines = tokensToStringset(space, "\r\n", content, contentlen);
#endif
  FREEMEMORY(space, content);
  *linecount=lines->noofstrings;
  csv=ALLOCMEMORY(space, NULL, stringset_t *, lines->noofstrings);

  for(i=0; i < lines->noofstrings; i++) {
    csv[i] = tokensToStringset(space, delim, lines->strings[i].str, lines->strings[i].len);
  }

  destructStringset(space, lines);	
  return csv;
}

void
writeY(char *filename, double *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%f\n", i, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUint(char *filename, Uint *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%d\n", i, Y[i]);
  }

  fclose(file);
  return;
}

void 
writeXYUint(char *filename, Uint *X, Uint *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%d\t%d\n", i, X[i], Y[i]);
  }

  fclose(file);
}


