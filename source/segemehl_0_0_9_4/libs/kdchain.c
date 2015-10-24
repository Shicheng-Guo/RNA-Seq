
/*
 *  kdchain.c
 *  implementation of kdchain
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 04/29/2008 07:01:30 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 85 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-18 15:34:44 +0100 (Tue, 18 Nov 2008) $
 *
 *  Id: $Id: kdchain.c 85 2008-11-18 14:34:44Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdchain.c $
 *  
 */

#include "manout.h"
#include "kdchain.h"
#include "mathematics.h"
#include "sufarray.h"
#include "container.h"
#include "kdseed.h"
#include "debug.h"
#include "karlin.h"
#include <assert.h>

int
comp_gpos(const void *a, const void *b) {
  gmatch_t *l;
  gmatch_t *r;

  l = (gmatch_t*) a;
  r = (gmatch_t*) b;

  if (l->p < r->p) return -1;
  if (l->p > r->p) return 1;

  return 0;
}

inline void
reportchaincoords(gmatch_t *dest,
                  gmatch_t *a,
                  gmatch_t *b,
                  int ovq,
                  int dmis,
                  int dins,
                  int ddel) {
      
 DBG("a: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(a), FEND_Q(a), 
      FSTART_S(a), FEND_S(a),
      a->scr, a->mat, a->mis, a->ins, a->del); 
      
  DBG("b: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(b), FEND_Q(b), 
      FSTART_S(b), FEND_S(b),
      b->scr, b->mat, b->mis, b->ins, b->del); 
  
  DBG("dest: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(dest), FEND_Q(dest), 
      FSTART_S(dest), FEND_S(dest),
      dest->scr, dest->mat, dest->mis, dest->ins, dest->del); 
      
  DBG ("ovl:%d dmis:%d, dins:%d, ddel:%d\n", ovq, dmis, dins, ddel);

  return;
}

void
initgchains(gchains_t *pi) {
    pi->size = 0;
    pi->bestscr = 0;
    pi->bestscrpos=0;
    pi->chains = NULL;
}

inline int
misfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (dds == ddq)  return dds;
  if (dds > ddq)   return ddq;
  if (ddq > dds)   return dds;

  return 0;
}


inline int
delfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (dds > ddq)  return dds - ddq;

  return 0;
}


inline int
insfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (ddq > dds)  return ddq - dds;

  return 0;
}

inline void
joinfrag(gmatch_t *dest, 
    gmatch_t *a, 
    gmatch_t *b) {

  int dmis=0,
      dins=0,
      ddel=0,
       ovq=0;
  

  if (a->scr > 0 && b->scr > 0) { 
    ovq  = abs(MIN(1, misfrag(a,b))-1);
    dmis = MAX(1, misfrag(a,b))-1;
    dins = insfrag(a,b);
    ddel = delfrag(a,b);
  } else {
    initMatch(dest); 
    return;
  }

  dest->scr = a->scr + b->scr - ovq - dmis - dins -ddel;
  dest->mat = a->mat + b->mat - ovq;
  dest->mis = a->mis + b->mis + dmis;
  dest->ins = a->ins + b->ins + dins;
  dest->del = a->del + b->del + ddel;

  dest->i   = a->i;
  dest->p   = a->p;
  
  if (dest->scr > 0 && FEND_Q(b) != FEND_Q(dest)) { 
    reportchaincoords(dest,a,b,ovq,dmis,dins,ddel);
  }
  
  if (dest->scr > 0 && FEND_S(b) != FEND_S(dest)) { 
    reportchaincoords(dest,a,b,ovq,dmis,dins,ddel);
  }

  dest->j   = FEND_Q(b);
  dest->q   = FEND_S(b);
  dest->checklen = a->checklen;
  dest->subject = a->subject;
  return;
}

inline void
reportfchain(void *space, 
    gchains_t *pi,  
    gmatch_t *e) {

  pi->chains = ALLOCMEMORY(space, pi->chains, gmatch_t, pi->size+1);
  memmove(&pi->chains[pi->size], e, sizeof(gmatch_t));
  if (e->scr > pi->bestscr) {
    pi->bestscr = e->scr;
    pi->bestscrpos = pi->size;
  }
  pi->size++;
  return;
}


gmatch_t*
greedyfchain( void *space,
    gmatch_t *F,
    Uint n,
    Uint *scr,
    Uint *pos,
    Uint *m) {


  int ovl;
  Uint i = 0;
  gmatch_t *h,
           *tmp1,
           *tmp2,
           *f,
           *e;
  gmatch_t *chains;
  gchains_t *pi=NULL;

  f = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  e = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  h = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  tmp1 = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  tmp2 = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  pi = ALLOCMEMORY(space, NULL, gchains_t, 1);

  initgchains(pi);
  initMatch(h);
  memmove(e, &F[0], sizeof(gmatch_t));

  for(i=0; i < n; i++)  {
    memmove(f, &F[i], sizeof(gmatch_t));
    ovl = FEND_S(e) - FSTART_S(f);

    if (FEND_Q(f) > FEND_Q(e) && FSTART_Q(f) < FSTART_Q(e)) {
      memmove(e, f, sizeof(gmatch_t));
    /*case 1: better start*/
    } else if (FEND_S(f) == FEND_S(e) && f->scr >= e->scr) {
      memmove(e, f, sizeof(gmatch_t));
    }
    /*case 2: positive overlap*/
    else if (FEND_S(f) > FEND_S(e) && ovl > 0) {
      joinfrag(tmp1, e, f);
      if(tmp1->scr > h->scr) {
        memmove(h,tmp1,sizeof(gmatch_t));
      } 
    }
    /*case 3: negative/zero overlap*/
    else if (FEND_S(f) > FEND_S(e) && ovl <= 0) {
      joinfrag(tmp1, h,f);
      joinfrag(tmp2, e,f);
      /*case 3a: bridge has better scr*/
      if (h->scr > 0 && FEND_S(f) > FEND_S(h) && tmp1->scr > tmp2->scr) {
        /*case 3aa: bridge of bridges*/
        if (FSTART_S(f) < FEND_S(h)) {
          memmove(e, h, sizeof(gmatch_t));
          memmove(h, tmp1, sizeof(gmatch_t));
        } else {
          memmove(e, tmp1, sizeof(gmatch_t));
          initMatch(h);
        }
      /*case 3b: accept d(e,f)*/
      } else if (tmp2->scr > e->scr) {
        memmove(e, tmp2, sizeof(gmatch_t));
        initMatch(h);
      /*case 3c: chain terminated*/
      } else {
        reportfchain(space, pi, e);
        memmove(e, f, sizeof(gmatch_t));
        initMatch(h);
      }
    } 
  }

  if(h->scr > e->scr) {
    reportfchain(space, pi, h);
  } else {
    reportfchain(space, pi, e);
  }

  FREEMEMORY(space, f);
  FREEMEMORY(space, e);
  FREEMEMORY(space, h);
  FREEMEMORY(space, tmp1);
  FREEMEMORY(space, tmp2);

  (*m) = pi->size;
  (*scr) = pi->bestscr;
  (*pos) = pi->bestscrpos;
  chains = pi->chains;

  FREEMEMORY(space, pi);
  
  return chains;
}

void
branch2match(Suffixarray *s,
    Container *C, 
    branch_t* b, 
    Uint noofbranches) {
  Uint i,
    k,
    idx;
  gmatch_t m;
  
  assert(noofbranches > 0);

  for(i=0; i < noofbranches; i++) {
    for(k=b[i].l; k <= b[i].r; k++) {
      m.p = s->suftab[k];
      m.scr = b[i].mat - b[i].mis - b[i].ins - b[i].del;
      m.mat = b[i].mat;
      m.del = b[i].del;
      m.ins = b[i].ins;
      m.mis = b[i].mis;
      m.i = i;
      m.j = i+b[i].mat+b[i].mis+b[i].ins-1;
      m.q = m.p+b[i].mat+b[i].mis+b[i].del-1;
      m.checklen = noofbranches;
      m.evalue = 0;
      m.scr = b[i].mat - b[i].mis - b[i].ins - b[i].del;
      idx = getMultiCharSeqIndex(s->seq, s->seq->sequences + m.p);
      m.subject = idx;
      if (b[i].mat > b[i].mis + b[i].ins + b[i].del) {
        bl_containerAdd(C, &m);
      }
    }
  }
  qsort(C->contspace, bl_containerSize(C), sizeof(gmatch_t), comp_gpos);
}

Container*
findfchains(void *space,
    Suffixarray *s,
    matchstem_t* M,
    Uint m,
    Uint t,
    unsigned char strict,
    int sigmatch,
    double lambda,
    double H,
    double K,
    double maxevalue) {

  Uint i,
    j,
    k,
    start,
    range,
    size,
    bestscr,
    bestscrpos,
    offset, 
    chainno,
    idx;

  Container P, 
              *C;
  gmatch_t p,
           *tmp,
           *ptr,
           *last,
           *G; 

  C = (Container *) malloc(sizeof(Container));
  bl_containerInit(C, 1000, sizeof(gmatch_t));
  bl_containerInit(&P, 1000, sizeof(gmatch_t));

  for(i=0; i < m; i++) {
    for(j=0; j < M[i].noofbranches; j++) {
      if(M[i].branches[j].r - M[i].branches[j].l < t) {
        for(k=M[i].branches[j].l; k <= M[i].branches[j].r; k++) {
          p.p = s->suftab[k];
          p.scr = M[i].branches[j].mat - M[i].branches[j].mis - M[i].branches[j].ins - M[i].branches[j].del;
          p.mat = M[i].branches[j].mat;
          p.del = M[i].branches[j].del;
          p.ins = M[i].branches[j].ins;
          p.mis = M[i].branches[j].mis;
          p.i   = i;
          p.j   = p.i+M[i].branches[j].mat+M[i].branches[j].mis+M[i].branches[j].ins-1;
          p.q   = p.p+M[i].branches[j].mat+M[i].branches[j].mis+M[i].branches[j].del-1;
          p.checklen = m;
          idx = getMultiCharSeqIndex(s->seq, s->seq->sequences + p.p);
          p.subject = idx;
          if(M[i].branches[j].mat > M[i].branches[j].mis + M[i].branches[j].ins + M[i].branches[j].del) {
            bl_containerAdd(&P, &p);
          }
        }
      }
    }
  }

  if(bl_containerSize(&P)) {
    qsort(P.contspace, bl_containerSize(&P), sizeof(gmatch_t), comp_gpos);

    ptr = (gmatch_t*)bl_containerGet(&P,0);
    last = ptr;
    start = 0;

    for(i=1; i < bl_containerSize(&P); i++) {
      ptr = (gmatch_t*) bl_containerGet(&P,i);

      if ((ptr->p - last->p) > last->scr+RINTERVALSIZE) {
        tmp = bl_containerGet(&P,start);
        range = (last->p - tmp->p)+tmp->scr+1;
        size = i-1-start; 
        DBGL(5, "interval start %d-%d[%d,%d]\n", range, size, start, i-1);
        offset = tmp->p;

        G = greedyfchain(space, tmp, size, &bestscr, &bestscrpos, &chainno);
        DBGL(5, "sigmatch: %d, bestscr:%d\n", sigmatch, G[bestscrpos].scr); 
        if (G[bestscrpos].scr > 0)
          G[bestscrpos].evalue = evalue(lambda, K, spacemult(m, s->numofsuffixes, H, K), G[bestscrpos].scr);

        if (strict) { 
          if ((G[bestscrpos].scr >= sigmatch 
                && G[bestscrpos].evalue < maxevalue)) { /*G[bestscr] >= 12*/
            bl_containerAdd(C, &G[bestscrpos]);
            DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
          }
        } else {
          if ((G[bestscrpos].scr >= sigmatch 
                || G[bestscrpos].evalue < maxevalue)) {
            bl_containerAdd(C, &G[bestscrpos]);
            DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
          }
        }

        bestscrpos = 0;
        bestscr = 0;
        start = i; 
        free(G);
      }

      DBGL(5, "pos:%d mat:%d mis:%d ins:%d del:%d\n", ptr->p, ptr->scr, 
          ptr->mis, ptr->ins, ptr->del); 
      last = ptr;
    }

    tmp = bl_containerGet(&P,start);
    range = (last->p - tmp->p)+tmp->scr+1;
    size = i-1-start; 
    DBGL(5, "interval start %d-%d[%d,%d]\n", range, size, start, i-1);
    offset = tmp->p;

    G = greedyfchain(space, tmp, size, &bestscr, &bestscrpos, &chainno);
    DBGL(5, "sigmatch: %d, bestscr:%d\n", sigmatch, G[bestscrpos].scr); 
    if (G[bestscrpos].scr > 0)
      G[bestscrpos].evalue = evalue(lambda, K, spacemult(m, s->numofsuffixes, H, K), G[bestscrpos].scr);

    if (strict) { 
      if (G[bestscrpos].scr > 0 && (G[bestscrpos].scr >= sigmatch 
            && G[bestscrpos].evalue < maxevalue)) { /*G[bestscr] >= 12*/
        bl_containerAdd(C, &G[bestscrpos]);
        DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
      }
    } else {
      if (G[bestscrpos].scr > 0 && (G[bestscrpos].scr >= sigmatch 
            || G[bestscrpos].evalue < maxevalue)) {
        bl_containerAdd(C, &G[bestscrpos]);
        DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
      }
    }


    bestscrpos = 0;
    bestscr = 0;
    start = i; 
    FREEMEMORY(space, G);
  }
  bl_containerDestruct(&P, NULL);
  return C;
}


