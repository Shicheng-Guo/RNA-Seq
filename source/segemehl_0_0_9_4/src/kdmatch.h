#ifndef KDMATCH_H
#define KDMATCH_H

/*
 *
 *	relaxed.h
 *  Declarations for kdmatch
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 12/26/07 02:41:35 CST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: kdmatch.h 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/src/kdmatch.h $
 *  
 */

#include "charsequence.h"
#include "kdseed.h"


#define SCORE(M, X, I, D)       (M) - (X) - (I) - (D)
#define MATCHES(S, X, I, D)     (S) + (X) + (I) + (D)
#define LEN_Q(S, X, I, D)       (S) + (X) + (I) + (D) + (X) + (I)  
#define LEN_S(S, X, I, D)       (S) + (X) + (I) + (D) + (X) + (D)


void
matchkdseed(void *space, 
    Suffixarray *sarray, 
    fasta_t *reads, 
    Uint minsize,
    char *outfile,
    Uint *counter,
    unsigned char silent,
    Uint s_ext,
    Uint p_mis,
    Uint Xoff,
    Uint k_p,
    Uint rep_type,
    Uint hitstrategy,
    Uint bedist,
    unsigned char showalignment,
    double maxevalue,
    int acc,
    Uint M,
    unsigned char matchingstat,
    FILE *dev,
    FILE *uninformativedev);

void kmismatch(void *space, 
    Suffixarray *s, 
    fasta_t *reads, 
    Uint k, 
    Uint *counter, 
    Uint rep_type,
    unsigned char silent,
    FILE *dev);

#endif
