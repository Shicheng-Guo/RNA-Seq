
/*
 *  matchfilesfields.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 01.03.2011 17:37:06 CET
 *  
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "alignment.h"
#include "debug.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "sort.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "zran.h"
#include "nw.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"


/*------------------------- bl_matchfileGetQname ------------------------------
 *    
 * @brief access query name of hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetQname(stringset_t *fields, unsigned char fmt) {

  if (fields->noofstrings < 4) return 0;

  switch(fmt) {
  case SAM:
    return fields->strings[0].str;
    break;

  default:
    DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);

  }

}



/*--------------------------- bl_matchfileGetRNext ---------------------------
 *    
 * @brief get mate chromosome
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetRNext ( stringset_t *fields, unsigned char fmt )
{
  
  if (fields->noofstrings < 6) return 0;

  switch(fmt) {
  case SAM:
    return fields->strings[6].str;
    break;

  default:
    DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);

  }
	return NULL;
}


/*--------------------------- bl_matchfileGetPNext ---------------------------
 *    
 * @brief get mate position
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetPNext (stringset_t *fields, unsigned char fmt)
{
  if (fields->noofstrings < 7) return 0;

  switch(fmt) {
  case SAM:
    return atoi(fields->strings[7].str);
    break;

  default:
    DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }
  return 0;
}

/*------------------------- bl_matchfileGetFlag -------------------------------
 *    
 * @brief access flag of hit
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetFlag(stringset_t *fields, unsigned char fmt) {

  if (fields->noofstrings < 4) return 0;

  switch(fmt) {
  case SAM:
    return atoi(fields->strings[1].str);
    break;

  default:
    DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);

  }

}

/*--------------------------- bl_matchfileIsHeader ---------------------------
 *    
 * @brief check header
 * @author Steve Hoffmann 
 *   
 */
 
unsigned char
bl_matchfileIsHeader (char *buffer, Uint len, unsigned char fmt)
{
  char *samhtags[] = {"@HD","@SQ","@RG", "@PG", "@CO"};
  int i, n = 5;

  switch(fmt) {
    case SAM:
      for(i=0; i < n; i++) {
   //     fprintf(stderr, "cmp tag '%s' with '%s'\n", samhtags[i], buffer);
        if(!strncmp(samhtags[i], buffer, 3)) return 1;
      }
      break;
    default:
      DBGEXIT("Unkown format (%d). Exit forced!\n", fmt);
  }
	
  return 0;
}

/*------------------------- bl_matchfileGetStartPos --------------------------
 *    
 * @brief access start position of hit
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetStartPos(stringset_t *fields, unsigned char fmt) {

  if(fields->noofstrings < 4) return 0;

  switch(fmt) {
    case SAM:
      return atoi(fields->strings[3].str);
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }
}

/*-------------------------- bl_matchfileGetEndPos ---------------------------
 *    
 * @brief access end position of hit
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetEndPos(stringset_t *fields, unsigned char fmt) {

  if(fields->noofstrings < 6) return 0;
  switch(fmt) {
    case SAM:
      return atoi(fields->strings[3].str) + 
        bl_cigarGetAlignLen(fields->strings[5].str) - 1;
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt); 
  }
}

/*--------------------------- bl_matchfileGetRead ----------------------------
 *    
 * @brief access read sequence of hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetRead(stringset_t *fields, unsigned char fmt) {

  if(fields->noofstrings < 10) return NULL;

  switch(fmt) {
    case SAM:
      return fields->strings[9].str;
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }
}

/*--------------------------- bl_matchfileGetQual ----------------------------
 *    
 * @brief access quality string for hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetQual(stringset_t *fields, unsigned char fmt) {
  Uint slen = 0;

  if(fields->noofstrings < 10) return NULL;
  switch(fmt) {
    case SAM:

      if (fields->strings[10].str[0] == '*' &&
          fields->strings[10].len == 1) {
        slen = strlen(fields->strings[9].str);
        fields->strings[10].str = 
          ALLOCMEMORY(space, fields->strings[10].str, char, slen+1);
        memset(fields->strings[10].str, 'b', slen);
        fields->strings[10].str[slen] = 0;
      }
      return fields->strings[10].str;
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }
}

/*-------------------------- bl_matchfileGetStrand ---------------------------
 *    
 * @brief access strandiness of hit
 * @author Steve Hoffmann 
 *   
 */
 
char
bl_matchfileGetStrand(stringset_t *fields, unsigned fmt) {
  char strands[] = {'+','-'};

  if(fields->noofstrings < 6) return -1;
  switch(fmt) {
    case SAM:
      return strands[(atoi(fields->strings[1].str) & 0x10) >> 4];
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }

  return 0;
}


/*--------------------------- bl_matchfileGetChrom ---------------------------
 *    
 * @brief access chromosome of hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetChrom(stringset_t *fields, unsigned fmt) {

  if(fields->noofstrings < 4) return 0;
  switch(fmt) {
    case SAM:
      return fields->strings[2].str;
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }

  return NULL;
}

/*------------------------- bl_matchfileGetMatchCnt --------------------------
 *    
 * @brief access number of parallel hits
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetMatchCnt(stringset_t *fields, unsigned char fmt) {
  Uint i, xn=1;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(fields->strings[i].len > 5 && 
            (strncmp(fields->strings[i].str, "XN:i:", 5) == 0 || 
             strncmp(fields->strings[i].str, "NH:i:", 5) == 0)) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*--------------------------- bl_matchfileGetEdist ---------------------------
 *    
 * @brief access edist of hit
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetEdist(stringset_t *fields, unsigned char fmt) {
  Uint i, nm = 0;

  switch(fmt) {
    case SAM:
     if(fields->noofstrings < 12) return 0;
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "NM:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          nm= atoi(&fields->strings[i].str[5]);
        }
      } 
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }

  return nm;
}

/*---------------------------- bl_matchfileGetAln ----------------------------
 *    
 * @brief access alignment string for hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetAln(stringset_t *fields, unsigned char fmt) {

  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
      return bl_cigarGetAlignString(fields->strings[5].str);
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }

  return NULL;
}

/*------------------------ bl_matchfileGetPrevPos -------------------------
 *    
 * @brief access the next split hit position
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetDiffString(stringset_t *fields, unsigned char fmt) {
  Uint i;
  char *res = NULL;

  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "MD:Z:", 5) == 0 &&
            fields->strings[i].len > 5) {
          res = bl_mdGetDiffString(&fields->strings[i].str[5]);
          return res;
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
  }

  return NULL;
}

/*------------------------- bl_matchfileGetNextPos --------------------------
 *    
 * @brief access the previous hit position
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetPrevPos(stringset_t *fields, unsigned char fmt) {
  Uint i, xn=1;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        
        if(strncmp(fields->strings[i].str, "XU:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*------------------------- bl_matchfileGetNextPos --------------------------
 *    
 * @brief access the previous hit position
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetNextPos(stringset_t *fields, unsigned char fmt) {
  Uint i, xn=1;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XV:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*------------------------- bl_matchfileGetPrevChr --------------------------
 *    
 * @brief access number of parallel hits
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetPrevChr(stringset_t *fields, unsigned char fmt) {
  Uint i;
  char *chr = NULL;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XP:Z:", 5) == 0 &&
            fields->strings[i].len > 5) {
          //chr = ALLOCMEMORY(space, NULL, char, fields->strings[i].len-4);
          //memmove(chr, &fields->strings[i].str[5], fields->strings[i].len-5);
          //chr[fields->strings[i].len-5] = 0;
          chr = &fields->strings[i].str[5];
        }
      }

      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return chr;
}

/*------------------------- bl_matchfileGetNextChr --------------------------
 *    
 * @brief access chromosome for next split hit
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGetNextChr(stringset_t *fields, unsigned char fmt) {
  Uint i;
  char *chr = NULL;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XC:Z:", 5) == 0 &&
            fields->strings[i].len > 5) {
          //chr = ALLOCMEMORY(space, NULL, char, fields->strings[i].len-4);
          //memmove(chr, &fields->strings[i].str[5], fields->strings[i].len-5);
          //chr[fields->strings[i].len-5] = 0;
          chr = &fields->strings[i].str[5];
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return chr;
}

/*------------------------ bl_matchfileGetSplitStart ------------------------
 *    
 * @brief access the start of the split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetSplitStart(stringset_t *fields, unsigned char fmt) {
  Uint i, xn =0;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XX:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*------------------------ bl_matchfileGetSplitEnd --------------------------
 *    
 * @brief access the end of split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetSplitEnd(stringset_t *fields, unsigned char fmt) {
  Uint i, xn = 0;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XY:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*----------------------- bl_matchfileGetSplitNumber ------------------------
 *    
 * @brief access the end of split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetSplitNumber(stringset_t *fields, unsigned char fmt) {
  Uint i, xn = 0;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XQ:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          xn = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return xn;
}

/*----------------------- bl_matchfileGetPrevFlag ------------------------
 *    
 * @brief get flags of the prev split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetPrevFlag(stringset_t *fields, unsigned char fmt) {
  Uint i, fl = 0;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XS:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          fl = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return fl;
}

/*----------------------- bl_matchfileGetNextFlag ------------------------
 *    
 * @brief get flags for the next split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetNextFlag(stringset_t *fields, unsigned char fmt) {
  Uint i, fl = 0;
      
  switch(fmt) {
    case SAM:
      if(fields->noofstrings < 12) return 0;
 
      for(i=11; i < fields->noofstrings; i++) {
        if(strncmp(fields->strings[i].str, "XT:i:", 5) == 0 &&
            fields->strings[i].len > 5) {
          fl = atoi(&fields->strings[i].str[5]);
        }
      }
      break;
    default:
      DBGEXIT("Unknown format (%d). Exit forced!\n", fmt);
    }
            
  return fl;
}


