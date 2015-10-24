
/*
 *  kdseed.c
 *  getting k-diff seeds
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 04/11/2008 10:41:11 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 91 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-07 16:44:27 +0100 (Sun, 07 Dec 2008) $
 *
 *  Id: $Id: kdseed.c 91 2008-12-07 15:44:27Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdseed.c $
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "memory.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include "biofiles.h"
#include "vtprogressbar.h"
#include "sort.h"
#include "bitArray.h"
#include "vqueue.h"
#include "vstack.h"
#include "container.h"
#include <pthread.h>
#include "kdseed.h"
#include "info.h"
#include "debug.h"
#include <assert.h>
/*#include "uedist.h"*/

unsigned int rekCounter=0;
unsigned int update=0;

inline int
kdscore(branch_t *b) {

  if (b->mat  < b->mis + b->ins + b->del) return 0;  
  return (int)b->mat - ((int)b->ins + (int)b->mis + (int)b->del);
}

inline void
matchstemAddBranch(void *space, 
    matchstem_t *a, 
    Uint mat, Uint q, 
    Uint mis, Uint ins, Uint del, 
    Uint l, Uint r, Uint u, Uint v) {
  
  a->branches = realloc(a->branches, sizeof(branch_t)*(a->noofbranches+1));
  a->noofbranches++;  
  a->branches[a->noofbranches-1].mat = mat;
  a->branches[a->noofbranches-1].q = q;
  a->branches[a->noofbranches-1].mis = mis;
  a->branches[a->noofbranches-1].ins = ins;
  a->branches[a->noofbranches-1].del = del;
  a->branches[a->noofbranches-1].l = l;
  a->branches[a->noofbranches-1].r = r;
  a->branches[a->noofbranches-1].u = u;
  a->branches[a->noofbranches-1].v = v;

  return;
}


inline void
matchstemModifyBranch(void *space, 
    matchstem_t *a, Uint k,
    Uint mat, Uint q, 
    Uint mis, Uint ins, Uint del, 
    Uint l, Uint r, Uint u, Uint v) {
 
  assert(a->noofbranches > k);

  a->branches[k].mat = mat;
  a->branches[k].q = q;
  a->branches[k].mis = mis;
  a->branches[k].ins = ins;
  a->branches[k].del = del;
  a->branches[k].l = l;
  a->branches[k].r = r;
  a->branches[k].u = u;
  a->branches[k].v = v;
  
  return;
}

void
matchstemDestruct(void *space, matchstem_t *M) {
  FREEMEMORY(space, M->branches);
}

/*-------------------------------- pushkdlcp ---------------------------------
 *    
 * @brief helper function to push lcps and singletons to a stack
 * @author Steve Hoffmann 
 *   
 */


inline void
pushkdlcp ( void *space, 
    VStack *vstack, 
    char *p, 
    Uint m, 
    Uint kp,
    kdiffm_t *data) {
  Uint i;
  bl_vstackPush(vstack, data);
  kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
  tmp->parent.a = tmp->child.a;
  tmp->parent.b = tmp->child.b;
  tmp->sptr++; tmp->kcnt++; tmp->del++;
  for(i = 1; i + data->kcnt <= kp && data->qptr + i < m; i++) {
    bl_vstackPush(vstack, data);
    kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
    tmp->parent.a = tmp->child.a;
    tmp->parent.b = tmp->child.b;
    tmp->qptr += i; tmp->kcnt += i; tmp->ins +=i;
  }

  return;
}

inline void
pushkdbranchesArr ( void *space,
    Suffixarray *s,
    VStack *vstack, 
    char *p, 
    Uint m,
    Uint kp,
    kdiffm_t *data, 
    PairUint pr) {

  Uint l,r,v,i;
  PairUint child;
  Lint *c;
  Uint count=0, lcp=0, j=0;

  child.a = 0;
  child.b = 0;

  c = getChildintervalsArr(space, s, pr.a, pr.b, &count);
  lcp = getlcpval(s, pr.a, pr.b);

  for(i=0; i < count; i++) {
    if(s->seq->sequences[ s->suftab[c[i*2]] + lcp] == p[data->qptr]){
      child.a = c[i*2];       
      child.b = c[i*2+1];
      break;
    }
  }
  
  for (v=0; v < count; v++) {
    l = c[v*2];
    r = c[v*2+1];
    if (l <= r) {
      if (l != child.a || r != child.b) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->qptr++; tmp->sptr++; tmp->kcnt++; tmp->mis++;
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = l; tmp->child.b = r;
      }
      if (data->qptr > 1) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->sptr++; tmp->kcnt++; tmp->del++; 
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = l; tmp->child.b = r;
      }
    }   
  }  

  /*Christian: i>0*/
  if (data->qptr > 1) {
    for(i = 1; i + data->kcnt <= kp && data->qptr + i < m; i++){
      
      for(j=0; j < count; j++) {
        if(s->seq->sequences[ s->suftab[c[j*2]] + lcp] == p[data->qptr+i]){
          child.a = c[j*2];       
          child.b = c[j*2+1];
          break;
        }
      }
      if (child.a <= child.b) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->qptr += i; tmp->kcnt += i; tmp->ins += i;
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = child.a; tmp->child.b = child.b;
      }
    }
  }

  free(c);
  return;
}


inline void
pushkdbranches ( void *space,
    Suffixarray *s,
    VStack *vstack, 
    char *p, 
    Uint m,
    Uint kp,
    kdiffm_t *data, 
    PairUint pr) {

  Uint l,r,v,i;
  Container *c;
  PairUint child;
  child = getCharInterval(space, s, pr.a, pr.b, 0, p[data->qptr]);
  c = getChildintervals(space, s, pr.a, pr.b);

  for (v=0; v < bl_containerSize(c); v++) {
    l = ((PairUint*)bl_containerGet(c,v))->a;
    r = ((PairUint*)bl_containerGet(c,v))->b;
    if (l <= r) {
      if (l != child.a || r != child.b) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->qptr++; tmp->sptr++; tmp->kcnt++; tmp->mis++;
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = l; tmp->child.b = r;
      }
      if (data->qptr > 1) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->sptr++; tmp->kcnt++; tmp->del++; 
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = l; tmp->child.b = r;
      }
    }   
  }  

  bl_containerDestruct(c, NULL);
  free(c);
  /*Christian: i>0*/
  if (data->qptr > 1) {
    for(i = 1; i + data->kcnt <= kp && data->qptr + i < m; i++){
      child = getCharInterval(space, s, pr.a, pr.b, 0, p[data->qptr+i]);      
      if (child.a <= child.b) {
        bl_vstackPush(vstack, data);
        kdiffm_t * tmp = (kdiffm_t *) bl_vstackTop(vstack);
        tmp->qptr += i; tmp->kcnt += i; tmp->ins += i;
        tmp->parent.a = pr.a; tmp->parent.b = pr.b;
        tmp->child.a = child.a; tmp->child.b = child.b;
      }
    }
  }
  return;
}


inline void 
kd_updatebranches(matchstem_t* b, kdiffm_t *data, unsigned int seedlen) {

  if(data->mis+data->ins+data->del == 0) {
    /*store the best branch*/
    b->branches = realloc(b->branches, sizeof(branch_t)*2);
    b->noofbranches = 2;
    matchstemModifyBranch(NULL, b, 1, data->mat, data->mat, data->mis, data->ins, data->del, 
        data->child.a, data->child.b, data->parent.a, data->parent.b);
    /*store the seed*/
    if(data->mat <= seedlen) {
      matchstemModifyBranch(NULL, b, 0, data->mat, data->mat, data->mis, data->ins, data->del, 
          data->child.a, data->child.b, data->parent.a, data->parent.b);
    }
  } else {  
    if(data->mat > data->ins+data->mis+data->del && data->mat-data->ins-data->mis-data->del > kdscore(&b->branches[1])) {
      b->branches = realloc(b->branches, sizeof(branch_t)*2);
      b->noofbranches = 2;
      matchstemModifyBranch(NULL, b, 1, data->mat, data->mat, data->mis, data->ins, data->del, 
          data->child.a, data->child.b, data->parent.a, data->parent.b);
    } else 
      if(data->mat > data->ins+data->mis+data->del && data->mat-data->ins-data->mis-data->del == kdscore(&b->branches[1])) {
      matchstemAddBranch(NULL, b, data->mat, data->mat, data->mis, data->ins, data->del, 
          data->child.a, data->child.b, data->parent.a, data->parent.b);
    }
  }
}



inline int
kd_matchlcp(void *space,
            Suffixarray *s, char *p, unsigned int m, 
            kdiffm_t *data, unsigned int cursufpos, 
            unsigned int lcplen, unsigned int seedlen, VStack *vstack,
            unsigned int maxdiff, int sext, int pmis, int xoff, matchstem_t *b) {

  unsigned int i;
  unsigned int fallbackmat=data->mat, fallbackmis=data->mis;
  unsigned int pnlty = 0;
  unsigned int seqlen = s->numofsuffixes;
  unsigned int lastmat=0, lastmis=0;
  char *cursuf;
  unsigned char  delimbreak=0;
  branch_t branch;

  /*to match inbetween the lcp*/
  cursuf = &s->seq->sequences[cursufpos];

  /*
   * we have to consider two indices:
   *    i iters through the lcp
   *    sptr shows the suffix pos (relative to cursuf) we start in the lcp
   *    qptr shows the pattern pos we start with
   *    we use "stacking" if alignment to the end of the query is
   *    not possible
   */

  for(i=0; i < lcplen 
      && cursufpos+data->sptr+i < seqlen 
      && data->qptr+i < m; i++) {
    
    if(cursuf[data->sptr+i] == s->seq->delim) break;

    if(data->kcnt < maxdiff) {
      data->qptr += i; data->sptr+=i; 
      pushkdlcp(space, vstack, p, m, maxdiff, data);
      data->qptr -= i; data->sptr-=i;
    }

    if(cursuf[data->sptr+i] == p[data->qptr+i]) {
      data->mat++;
      lastmat = data->mat;
      lastmis = data->mis;
      /* 
       * we reduce the overall penalty by the number of extensions
       * and keep track of the last pnlty free hit
       */
      pnlty = MAX(0, (int)(pnlty-sext));

      if(!pnlty) {
        fallbackmat = data->mat;
        fallbackmis = data->mis;
      }
      if (data->mis+data->ins+data->del == 0)
        kd_updatebranches(b, data, seedlen);
    } else {
      /*
       * check if mismatch violates the xoff constraint
       */
      data->mis++;

      if((pnlty+=pmis) > xoff) {
        break;
      }
    }
  }
  /*
   * we need a clean suffix if lcp extension failed
   * and we accept the move otherwise
   * */

  if(cursuf[data->sptr+i] == s->seq->delim) {
    i -= (i > 0) ? 1 : 0;
    delimbreak =1;
  }

  if(pnlty > xoff) {
    data->mat = fallbackmat;
    data->mis = fallbackmis;
  } else {
    data->qptr+=i;
    data->sptr+=i;
    if (data->mis+data->ins+data->del == 0) {
      kd_updatebranches(b, data, seedlen);
    } else {
      branch.mat = lastmat;
      branch.mis = lastmis;
      branch.ins = data->ins;
      branch.del = data->del;
      if (kdscore(&branch) > kdscore(&b->branches[1])) {
        b->branches = realloc(b->branches, sizeof(branch_t)*2);
        b->noofbranches=2;
        matchstemModifyBranch(space, b, 1, lastmat, lastmat, lastmis, data->ins, data->del, 
            data->child.a, data->child.b, data->parent.a, data->parent.b);
      } else {
        if (kdscore(&branch) == kdscore(&b->branches[1])){
          matchstemAddBranch(space, b, lastmat, lastmat, lastmis, data->ins, data->del, 
              data->child.a, data->child.b, data->parent.a, data->parent.b); 
        }
      }
    }
    if (delimbreak) {
      return 0;
    } else return 1;
  }

  if (delimbreak) {
    return 0;
  }

  kd_updatebranches(b, data, seedlen);

  return (!(pnlty>xoff));
}



inline matchstem_t*
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
    Uint qptr ) {

  Uint lcplen = 0,
       cursufpos = 0,
       seedlen=10000;
  VStack *vstack;
  PairUint birth;
  matchstem_t *b;
  kdiffm_t data, *tmp;

  data.kcnt = 0; data.mis = 0; data.mat = 0; data.ins = 0; data.del = 0;
  data.sptr = sptr; data.qptr = qptr;
  data.child.a = 1; data.child.b = 0;

  b = (matchstem_t *) malloc(sizeof(matchstem_t));
  b->noofbranches = 0;
  b->branches = malloc(sizeof(branch_t)*2);

  vstack = (VStack *) malloc(sizeof(VStack));
  bl_vstackInit(vstack, 1000, sizeof(kdiffm_t));

  data.l = getlcpval(s, iprime, jprime);
  data.parent.a = iprime; data.parent.b = jprime;

  if (data.l == data.sptr && data.l < m && iprime != jprime){
    data.child = getCharIntervalArr(space, s, iprime, jprime, 0, p[data.qptr]);

    if (data.ins+data.del+data.mis < maxdiff){
      pushkdbranches(space, s, vstack, p, m, maxdiff, &data, data.parent);
    }
  } 

  if (data.child.a > data.child.b){
    data.child.a = iprime;
    data.child.b = jprime;
  } else {
    data.mat++; data.sptr++; data.qptr++;
  }


  matchstemAddBranch(space, b, data.mat, data.mat, data.mis, data.ins, data.del, 
       data.child.a, data.child.b, data.parent.a, data.parent.b); 
  matchstemAddBranch(space, b, data.mat, data.mat, data.mis, data.ins, data.del, 
       data.child.a, data.child.b, data.parent.a, data.parent.b); 


  while(1) {
    while(data.qptr < m && data.child.a <= data.child.b) {
      data.l = getlcpval(s, data.child.a, data.child.b);
      cursufpos = s->suftab[data.child.a];
      /*
       * intermediate longest common prefixes (l > sptr)
       * singletons or (implicitly: lcps exceeding m)
       */
      if(data.l > data.sptr || data.l == 0) {
        /* if the lcp exceeds m we align and break*/
        if(data.l >= m + data.del - data.ins || data.l == 0) {
          //kd_align(s, p, m, cursufpos, &data, b, seedlen);
          kd_matchlcp(space, s, p, m, &data, cursufpos, m-data.del, seedlen, vstack,
                maxdiff, sext, pmis, xoff, b);
          /*full match triggers return otherwise 
           *evaluation of branches is continued*/
          if(data.mat==m) {
            bl_vstackDestruct(vstack, NULL);
            free(vstack);
            return b;
          } else { 
            break;
          }
        } else {
          lcplen = data.l - data.sptr;
          if(!kd_matchlcp(space, s, p, m, &data, cursufpos, lcplen, seedlen, vstack,
                maxdiff, sext, pmis, xoff, b) ){
            break;
          } 
        }
      } else {

      birth = getCharIntervalArr(space, s, data.child.a, 
          data.child.b, 0, p[data.qptr]);
      if(data.kcnt < maxdiff) {
        pushkdbranches(space, s, vstack, p, m, maxdiff, &data, data.child);
      }
      /*this character was a dead end*/
      if(birth.a > birth.b) {
        break;
      } 
      /*otherwise accept move*/
      data.parent.a = data.child.a;
      data.parent.b = data.child.b;
      data.child.a = birth.a;
      data.child.b = birth.b;
      data.sptr++;
      data.qptr++;
      data.mat++;
      kd_updatebranches(b, &data, seedlen);
    }
  }
  if(bl_vstackIsEmpty(vstack)) break;
  tmp = (kdiffm_t*) bl_vstackPop(vstack, NULL);
  memcpy(&data, tmp, sizeof(kdiffm_t));
  free(tmp);
}
bl_vstackDestruct(vstack, NULL);
free(vstack);
return(b);
}




inline matchstem_t*
kdseeds ( void *space,
    Suffixarray *s,
    char *p,
    Uint m,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint kp ) {
  Uint iprime,
  jprime,
  kprime,
  lprime,
  l,
  j,
  r,
  u,
  v,
  i, 
  c=0,
  q,
  ll;
  PairUint suflink;
  matchstem_t *b,
              *a;

  iprime = 0;
  jprime = s->numofsuffixes-1;


  b = kd_match(space, s, p, m, sext, pmis, xoff, kp, iprime, jprime, iprime, jprime, 0, 0); 
  a = ALLOCMEMORY(space, NULL, matchstem_t, m);

  a[0].branches = ALLOCMEMORY(space, NULL, branch_t, b->noofbranches-1); 
  memmove(a[0].branches, &b->branches[1], sizeof(branch_t)*(b->noofbranches-1));
  a[0].noofbranches = b->noofbranches-1;

  q = b->branches[0].q; 
  l = b->branches[0].l; r = b->branches[0].r;
  u = b->branches[0].u; v = b->branches[0].v;

  FREEMEMORY(space, b->branches);
  FREEMEMORY(space, b);  

  for(i=1; i < m; i++) { 

    ll = getlcpval(s, l, r);

    if(q == ll && q > 1) {
      suflink = getSuflink(s, l, r);
      kprime = l;
      lprime = r;
      iprime = suflink.a;
      jprime = suflink.b;
      c = q-1;
    } else {
      if(q <= 1) {
        kprime = 0;
        iprime = 0;
        lprime = s->numofsuffixes-1; 
        jprime = s->numofsuffixes-1;
        c = 0;
      } else {      
        kprime = u; 
        lprime = v;
        suflink = getSuflink(s, u, v);
        iprime = suflink.a;
        jprime = suflink.b;
        c = getlcpval(s, iprime, jprime);
      }
    }

    if (i < m) {
      b = kd_match(space, s, &p[i], m-i, sext, pmis, xoff, kp, 
          iprime, jprime, kprime, lprime, c, c);
 
      for(j=0; j < b->noofbranches; j++) {
        b->branches[j].mat += c;
      }
      
      a[i].branches = ALLOCMEMORY(space, NULL, branch_t, b->noofbranches-1); 
      memmove(a[i].branches, &b->branches[1], sizeof(branch_t)*(b->noofbranches-1));
      a[i].noofbranches = b->noofbranches-1;      
      
      q = b->branches[0].q + c; 
      l = b->branches[0].l; r = b->branches[0].r;
      u = b->branches[0].u; v = b->branches[0].v;

      FREEMEMORY(space, b->branches);
      FREEMEMORY(space, b);

    } else {
      break;
    }
  }

  return a;
}


inline void
dumpkdseeds(Suffixarray *s, matchstem_t *M, Uint m, char strand, Uint T) {
  Uint i,j,k;

  //MSG("kdseeds:\n");
  for(i=0; i < m; i++) {
    for(k=0; k < M[i].noofbranches; k++) {
      printf("%d %c ", M[i].branches[k].mat, strand);
      if(M[i].branches[k].mat > 0 && M[i].branches[k].r-M[i].branches[k].l <= T) {
        for(j=M[i].branches[k].l; j <= M[i].branches[k].r; j++) {
          printf("%d ", s->suftab[j]);
        }
      }
    }
    printf("\n");
    //printf("%d:%d(i:%d, d:%d, m:%d)-(%d..%d)\t%c\n", i, M[i].mat, M[i].ins, M[i].del, M[i].mis, M[i].l, M[i].r, strand);
    // printf("%d:%d-(%d..%d)\t", i, M[i].mat, s->suftab[M[i].l], s->suftab[M[i].r]);
  }
}

void kdcompare(TripleSint *a, branch_t *b, Uint m) {
  Uint i;

  for(i=0; i < m; i++) {
    if (b[i].mat != a[i].c || b[i].l != a[i].a || b[i].r != a[i].b) {
      NFO("failure at %d of %d", i, m);
    }
  }
}

/*----------------------------------- kmis ------------------------------------
 *    
 * @brief enumerates all matches in the suffix array with $k$ mismatches for a 
 *        given pattern $P$ of length $m$.
 * @author Steve Hoffmann 
 *   
 */

inline branch_t*
kmis (void *space,
    Suffixarray *s,
    char *P,
    Uint m,
    Uint k,
    Uint *noofmatches) {

  branch_t *matches=NULL;
  VQueue vqueue;
  Uint   i,
         matchno=0;
  char *cursuf;
  PairUint child;
  kmis_t data;
  Container *c;
  int min, lcp=0, llcp=0;

  data.l = 0; data.r = s->numofsuffixes-1;
  child.a = data.l;
  child.b = data.r;
  data.u = 0; data.v = s->numofsuffixes-1;
  data.p = 0; data.mis = 0;

  bl_vqueueInit(&vqueue, 1000, sizeof(kmis_t));
  bl_vqueueEnqueue(&vqueue, &data);

  while (!bl_vqueueIsEmpty(&vqueue)) {
    kmis_t *tmp = bl_vqueueDequeue(&vqueue, NULL);
    memcpy(&data, tmp, sizeof(kmis_t));
    free(tmp);
    llcp = data.p;

    /* if not singleton scan the string and enqueue alternatives */
    if (data.l < data.r) {
      while(1){

        lcp = getlcpval(s, data.l, data.r);
        min = MIN(lcp, (m-1));
        if (lcp > llcp+1) {
          for(i = llcp; i < lcp && i < m; i++) {
            cursuf = &s->seq->sequences[s->suftab[data.l]];
            if(P[i] != cursuf[i]) {
              data.mis++;
            }
          }
        }
        if (lcp > m-1) {
	  child.a = data.l;
	  child.b = data.r;
	  break;
	}
        c = getChildintervals(space, s, data.l, data.r);
        child = getCharInterval(space, s, data.l, data.r, 0, P[lcp]);

        data.u = data.l;
        data.v = data.r;

        for (i = 0; i < bl_containerSize(c); i++) {

          data.l = ((PairUint*)bl_containerGet(c,i))->a;
          data.r = ((PairUint*)bl_containerGet(c,i))->b;

          if (data.l <= data.r) {
            if ((data.l != child.a || data.r != child.b) &&
                data.mis + 1 <= k) {
              bl_vqueueEnqueue(&vqueue, &data);
              kmis_t *tmp = bl_vqueueFrontN(&vqueue, bl_vqueueSize(&vqueue)-1);
              tmp->p = lcp + 1;
              tmp->mis++;
            }
          }
        }       
        bl_containerDestruct(c, NULL);
        free(c);

        if (child.a >= child.b) {
          break;
        }

        llcp = lcp;
        data.l = child.a;
        data.r = child.b;
      } 
    } else {
      child.a = data.l;
      child.b = data.r;
      lcp = data.p;
    }

    if (child.a == child.b) {
      for(i = lcp; i < m; i++) {
        cursuf = &s->seq->sequences[s->suftab[child.a]];
        if(i+s->suftab[child.a] > s->numofsuffixes || P[i] != cursuf[i]) {
          data.mis++;
        }
        if (data.mis > k) break;
      }
      data.r = child.a;
      data.l = child.a;
    }
    if(data.mis <= k && child.a <= child.b) {
      cursuf = &s->seq->sequences[s->suftab[child.a]];
      matches = (branch_t *) realloc(matches, sizeof(branch_t) * (matchno+1));
      KMSTOREBRANCH(matches, matchno);
      matchno++;
    }
  }
  bl_vqueueDestruct(&vqueue, NULL);

  *noofmatches = matchno;
  return matches;
}



