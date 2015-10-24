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

char* readfile(void *, char *, Uint*);
stringset_t **readcsv(void *, char *, char*, Uint *);
void writeY(char *, double  *, Uint);
void writeXYUint(char *filename, Uint *X, Uint *Y, Uint len);


#endif
