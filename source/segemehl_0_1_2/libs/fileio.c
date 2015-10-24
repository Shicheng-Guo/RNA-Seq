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
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "stringutils.h"
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



void
bl_writeFileHeader(char *filename, char *header) {
  char *tmpfilename, *buffer;
  FILE *out, *in;
  size_t buffersize= 1024, len;

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  tmpfilename = bl_getTempFile("headerwrite",11);
  
  out = fopen(tmpfilename, "w");
  
  if(!out) {
    fprintf(stderr, "Couldnt open file %s for writing. Exit forced.", tmpfilename);
    exit(-1);
  }

  fprintf(out, "%s\n", header);
  fclose(out);

  out = fopen(tmpfilename, "a");
  in = fopen(filename, "r");
  
  if(!in) {
    fprintf(stderr, "Couldnt open file %s for reading. Exit forced.", filename);
    exit(-1);
  }

  while((len = fread(buffer, 1, buffersize, in)) > 0) {
    fwrite(buffer, 1, len, out);
  }

  fclose(out);
  fclose(in);
  rename(tmpfilename, filename);

  return;
}


char* 
bl_replacenonalphanum(char *s, Uint len) {
  Uint i, u=0, lastalphanum=0;
  int ch;
  char *new = ALLOCMEMORY(NULL, NULL, char, len+1);

  for(i=0; i < len; i++) {
    ch = s[i];
    if(((Uint)((ch | 0x20) - 'a') < 26u)|| ((Uint)(ch-'0') < 10u)) {
      lastalphanum = u;
      new[u++] = ch;
    } else {
      if (u > 0)
      new[u++] = '_';
    }
  }

  new[lastalphanum+1] ='\0';
  return new;
}



/*------------------------------ bl_getTempFile ------------------------------
 *    
 * @brief get a temporary file
 * @author Steve Hoffmann 
 *   
 */

char *
bl_getTempFile(char *tmp, Uint tmplen)
{

  int res;
  char *fname=NULL;

  fname = ALLOCMEMORY(NULL, NULL, char, tmplen+11);
  
  if(tmplen > 0)
    sprintf(fname, "%sXXXXXX", tmp);
  else
    sprintf(fname,"XXXXXX");

  if ((res = mkstemp(fname)) == -1) {
    fprintf(stderr, "Error in creating temporary file '%s'. Exit forced.\n",
            fname);
    exit(-1);
  }
  if (close(res) == -1){
    fprintf(stderr, "Error in closing temporary file '%s'. Exit forced.\n",
            fname);
    exit(-1);   
  }
  return fname;
}

int
bl_UnixSortMerge(void *space, char **filenames, Uint nooffiles, 
    const char *fieldstring, char *outfile) {
  int ret;
  Uint i, filenamestringpos;
  char *prg = "LC_COLLATE=C sort";
  char *cmd;
  char *filenamestring = NULL;
  char *tempfilename;

  tempfilename = bl_getTempFile("sort", 4);

  filenamestringpos = 0;
  for(i = 0; i < nooffiles; i++) {
   
    filenamestring = 
      ALLOCMEMORY(space, filenamestring, char, filenamestringpos+strlen(filenames[i])+2);
    
    memmove(&filenamestring[filenamestringpos], filenames[i], strlen(filenames[i]));
    filenamestringpos += strlen(filenames[i]);
    filenamestring[filenamestringpos] = ' ';
    filenamestringpos++;
    filenamestring[filenamestringpos] = 0;
  }


  cmd = ALLOCMEMORY(space, NULL, char, strlen(prg) + strlen(fieldstring) 
      + strlen(filenamestring) + strlen(outfile) + 6);
  sprintf(cmd, "%s -m %s %s > %s", prg, fieldstring, filenamestring, outfile);
  ret = system(cmd);

  return ret;
}

int
bl_rm(void *space, char *filename) {
  int ret=0;
  char *prg = "rm";
  char *cmd;

  cmd = calloc(strlen(prg) + strlen(filename) + 10, 1);
 
  sprintf(cmd, "%s -f %s", prg, filename);
  system(cmd);
 
  free(cmd);
  return ret;
}

int
bl_UnixSort(void *space, char *filename, const char *fieldstring) {
  int ret=0;
  char *prg = "LC_COLLATE=C sort";
  char *cmd;
  char *tempfilename;

  tempfilename = bl_getTempFile("sort", 4);

  cmd = calloc(strlen(prg) + strlen(fieldstring) 
      + strlen(filename) + strlen(tempfilename) + 10, 1);
 
  sprintf(cmd, "%s -o %s %s %s", prg, tempfilename, fieldstring, filename);
  system(cmd);
  
  bl_rm(space, filename);
  rename(tempfilename, filename);

  free(cmd);
  free(tempfilename);
  return ret;
}

char* dirname(const char *filename) {
    char *s;
    s=strrchr(filename, (int)'/');
    if(s && *s)
      *s = '\0';

    return s;
}

int
bl_fileprefixlen(char *filename) {
  Uint i, suf;
  
  suf = strlen(filename);
  for(i=1; i < strlen(filename); i++) {
    if(filename[i] == '.') {
      suf = i;
    }
  }
  return suf;
}

char *
bl_basename (const char *name)
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


void
bl_fnreplace(char *filename, char oldchar, char newchar, Uint nreplace) {
  int ch, n=0;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    if(ch == oldchar) {
      fseek(fp, -1, SEEK_CUR);
      fputc(newchar, fp);
      n++;
    }
    if(nreplace == n) break;
  }


  fclose(fp);
  return;
}

void
bl_freplacearr(char *filename, char* oldchars, char* newchars, Uint len, char stop) {
  int ch, i;
  char oldchar;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    for(i=0; i < len; i++) {
      oldchar = oldchars[i];
      if(ch == oldchar) {
        fseek(fp, -1, SEEK_CUR);
        fputc(newchars[i], fp);
        break;
      }
    }

    if(ch == stop) {
      fseek(fp, -1, SEEK_CUR);
      fputc(' ', fp);
      break;
    }
  }

  fclose(fp);
  return;
}

void
bl_freplace(char *filename, char oldchar, char newchar, char stop) {
  int ch;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    if(ch == oldchar) {
      fseek(fp, -1, SEEK_CUR);
      fputc(newchar, fp);
    }
    if(ch == stop) {
      fseek(fp, -1, SEEK_CUR);
      fputc('\n', fp);
      break;
    }
  }

  fclose(fp);
  return;
}

void
bl_freplacestr(char *filename, char *str, Uint len, char stop){
  int i = 0;
  char ch;
  FILE *fp;

  fp = fopen(filename, "rb+");  
  if (!fp) {
    fprintf(stderr, "Couldn't open file '%s'. Exit forced.\n", filename);
    exit(EXIT_FAILURE);
  }

  while((ch = fgetc(fp)) != EOF){
    if (ch == stop){
      break;
    }
    fseek(fp, -1, SEEK_CUR);
    fputc(str[i%len], fp);
    i++;
  }

  fclose(fp);
  return;
}

int
bl_fgets(void *space, FILE *fp, char **str) {
  char ch, *buffer;
  Uint buffersize = 100;
  Uint len = 0;
  
  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  
  while((ch=getc(fp)) != EOF && ch != '\n') {
    if(len == buffersize - 1) {
      buffersize = 2 * buffersize + 1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }
    buffer[len++] = (char) ch;
  }

  if(ch == EOF) return EOF;

  buffer[len] = '\0';
  *str = buffer;

  return len;
}


char* 
readfile(void* space, char* filename, Uint* strlen) {

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
writeY(char *filename, double *Y, Uint len, Uint xoff, Uint yoff) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=yoff; i < len; i++) {
    fprintf(file,"%d\t%f\n", i+xoff, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUint(char *filename, Uint *Y, Uint len, Uint xoff, Uint yoff) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=yoff; i < len; i++) {
    fprintf(file,"%d\t%d\n", i+xoff, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUintNorm(char *filename, Uint *Y, Uint len, Uint off) {
  FILE *file;
  Uint i, norm=0;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }
  for(i=0; i < len; i++) {
    norm += Y[i];
  }

  for(i=off; i < len; i++) {
    fprintf(file,"%d\t%f\n", i, (double)Y[i]/norm);
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

void
writeXYZ(char *filename, double *X, double *Y, double *Z, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%f\t%f\t%f\n", X[i], Y[i], Z[i]);
  }

  fclose(file);
}

