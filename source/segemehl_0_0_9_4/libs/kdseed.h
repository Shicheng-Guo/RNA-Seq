#ifndef KDSEEDS_H
#define KDSEEDS_H

/*
 *
 *	kdseed.h
 *  gettin k-diff seeds
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 04/13/2008 12:05:48 AM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 77 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-17 13:16:59 +0100 (Mon, 17 Nov 2008) $
 *
 *  Id: $Id: kdseed.h 77 2008-11-17 12:16:59Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdseed.h $
 */

#define MSTEMSTOREBRANCH(B, P)\
                    B[P].mat = data.mat;\
                    B[P].q = data.mat;\
                    B[P].mis = data.mis;\
                    B[P].ins = data.ins;\
                    B[P].del = data.del;\
                    B[P].l = data.child.a;\
                    B[P].r = data.child.b;\
                    B[P].u = data.parent.a;\
                    B[P].v = data.parent.b

#define MSTEMSTOREBRANCHP(B, P)\
                    B[P].mat = data->mat;\
                    B[P].q = data->mat;\
                    B[P].mis = data->mis;\
                    B[P].ins = data->ins;\
                    B[P].del = data->del;\
                    B[P].l = data->child.a;\
                    B[P].r = data->child.b;\
                    B[P].u = data->parent.a;\
                    B[P].v = data->parent.b

#define KDSTOREBRANCHP(B, P)\
                     B[P].mat = data->mat;\
                    B[P].q = data->mat;\
                    B[P].mis = data->mis;\
                    B[P].ins = data->ins;\
                    B[P].del = data->del;\
                    B[P].l = data->child.a;\
                    B[P].r = data->child.b;\
                    B[P].u = data->parent.a;\
                    B[P].v = data->parent.b

#define KDSTOREBRANCHEXTP(B, P)\
                    B[P].mat = lastmat;\
                    B[P].q = lastmat;\
                    B[P].mis = lastmis;\
                    B[P].ins = data->ins;\
                    B[P].del = data->del;\
                    B[P].l = data->child.a;\
                    B[P].r = data->child.b;\
                    B[P].u = data->parent.a;\
                    B[P].v = data->parent.b


#define KDSTOREBRANCH(B, P)\
                     B[P].mat = data.mat;\
                    B[P].q = data.mat;\
                    B[P].mis = data.mis;\
                    B[P].ins = data.ins;\
                    B[P].del = data.del;\
                    B[P].l = data.child.a;\
                    B[P].r = data.child.b;\
                    B[P].u = data.parent.a;\
                    B[P].v = data.parent.b

#define KDSTOREBRANCHEXT(B, P)\
                    B[P].mat = lastmat;\
                    B[P].q = lastmat;\
                    B[P].mis = lastmis;\
                    B[P].ins = data.ins;\
                    B[P].del = data.del;\
                    B[P].l = data.child.a;\
                    B[P].r = data.child.b;\
                    B[P].u = data.parent.a;\
                    B[P].v = data.parent.b

#define KDSTOREBRANCHEXTP(B, P)\
                    B[P].mat = lastmat;\
                    B[P].q = lastmat;\
                    B[P].mis = lastmis;\
                    B[P].ins = data->ins;\
                    B[P].del = data->del;\
                    B[P].l = data->child.a;\
                    B[P].r = data->child.b;\
                    B[P].u = data->parent.a;\
                    B[P].v = data->parent.b

#define KMSTOREBRANCH(B,P)\
                    B[P].mat = m-data.mis;\
                    B[P].q   = m-data.mis;\
                    B[P].mis = data.mis;\
                    B[P].ins = 0;\
                    B[P].del = 0;\
                    B[P].l = data.l;\
                    B[P].r = data.r;\
                    B[P].u = data.u;\
                    B[P].v = data.v

#define KDSCORE(B, P)\
                    kdscore(B[P].mat,B[P].mis,B[P].ins,B[P].del)

#include "basic-types.h"
#include "vstack.h"
#include "sufarray.h"



typedef struct branch_s {
    Uint mis;
    Uint mat;
    Uint q;
    Uint ins;
    Uint del;
    Uint l;
    Uint r;
    Uint u;
    Uint v;
} branch_t;

typedef struct matchstem_s {
  
  branch_t *branches;
  Uint noofbranches;
} matchstem_t;



typedef struct {
  Uint sptr;
  Uint qptr;
  Uint l;
  Uint kcnt;
  Uint mat;
  Uint mis;
  Uint ins;
  Uint del;
  PairUint child;
  PairUint parent;
} kdiffm_t;

typedef struct {
  Uint l;
  Uint r;
  Uint u;
  Uint v;
  Uint p;
  Uint mis;
} kmis_t;


extern  int
kdscore(branch_t *b);

extern  void
dumpkdseeds(Suffixarray *s, matchstem_t *M, Uint m, char strand, Uint T);

matchstem_t*
kdseeds ( void *space,
          Suffixarray *s,
          char *p,
          Uint m,
          Uint sext,
          Uint pmis,
          Uint xoff,
          Uint kp );
   

extern  branch_t* 
kdiff ( void *space, 
	Suffixarray *s,
	char *p,
	Uint m,
	Uint kp,
	Uint iprime,
	Uint jprime,
	Uint kprime,
	Uint lprime,
	Uint sptr,
	Uint qptr );


extern int
kd_matchlcp(void *space,
            Suffixarray *s, char *p, unsigned int m, 
            kdiffm_t *data, unsigned int cursufpos, 
            unsigned int lcplen, unsigned int seedlen, VStack *vstack,
            unsigned int maxdiff, int sext, int pmis, int xoff, matchstem_t *b);
extern void
pushkdbranches ( void *space,
		 Suffixarray *s,
		 VStack *vstack, 
		 char *p, 
		 Uint m,
		 Uint kp,
		 kdiffm_t *data, 
		 PairUint pr);

extern void
pushkdlcp ( void *space, 
	    VStack *vstack, 
	    char *p, 
	    Uint m, 
	    Uint kp,
	    kdiffm_t *data);

extern matchstem_t*
kdseeds ( void *space,
    Suffixarray *s,
    char *p,
    Uint m,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint kp );

extern matchstem_t*
kd_match ( void *space, 
    Suffixarray *s,
    char *p,
    Uint m,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint maxdiff,
    Uint iprime,
    Uint jprime,
    Uint kprime,
    Uint lprime,
    Uint sptr,
    Uint qptr );

extern void 
kd_updatebranches(matchstem_t* b, kdiffm_t *data, unsigned int seedlen);

void matchstemDestruct(void *space, matchstem_t *M);

  branch_t* kmis (void *space, Suffixarray *s, char *P, Uint m, Uint k, Uint *noofmatches);
void kdcompare(TripleSint *a, branch_t *b, Uint m);

#endif
