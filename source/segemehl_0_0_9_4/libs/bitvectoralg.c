
/*
 *  bitvectoralg.c
 *  implementation of Gene Myers 
 *  bitvector algorithm
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 05/23/2008 06:12:54 PM CEST
 *  
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "bitVector.h"
#include "memory.h"
#include "alignment.h"


/*-------------------------------- alphacheck ---------------------------------
 *    
 * @brief alphabet check
 * @author Steve Hoffmann 
 *   
 */

unsigned char 
alphacheck (char c) {
  if (   c == 'A' || c == 'a' || c == 'C' || c == 'c' 
      || c == 'T' || c == 't' || c == 'G' || c == 'g') {
    return 1;
  }
  return 0;
}


/*-------------------------------- encodetab ---------------------------------
 *    
 * @brief brutal ascii encoding
 * @author Steve Hoffmann 
 *   
 */

Uint*
encodetab(char *alphabet, Uint asize) {
    Uint i;
    Uint *tab;

    tab = ALLOCMEMORY(space, NULL, Uint, 255);
    memset(tab, 0, sizeof(Uint)*255);
    for(i=0; i < asize; i++) {
        tab[(Uint)alphabet[i]] = i;
    }
    return tab;
} 


/*---------------------------------- getpeq ----------------------------------
 *    
 * @brief returns pattern mask for each char in alphabet
 * @author Steve Hoffmann 
 *   
 */
 
bitvector*
getpeq(void *space,
    char *query, 
    Uint qlen,
    char *alphabet,
    Uint asize,
    Uint *enctab) {

  bitvector* peq;
  Uint i,j, wordno;

  wordno = qlen/BITVECTOR_WORDSIZE;
  wordno++;
  peq = ALLOCMEMORY(space, NULL, bitvector*, asize);

  for(i=0; i < asize; i++) {
    peq[i] = initbitvector(space, BITVECTOR_WORDSIZE*wordno);
    setbitvector(peq[i], BITVECTOR_WORDSIZE*wordno, 0);
    for(j=0; j < qlen; j++) {
      if (alphacheck(query[j]) && enctab[(Uint)query[j]] == i) {
        bitvector_setbit(peq[i], j, 1);
      } 
    }
  }
  return peq;
}


/*------------------------------ myersbitvector ------------------------------
 *    
 * @brief approx. string matching: calculate min edist of query and subject
 * @author Steve Hoffmann 
 *   
 */
 
PairSint
myersbitvector(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq) {

  bitvector 
  Pv,
  Mv,
  Eq;
  PairSint res;

  Uint score=qlen,
       i,
       j,
       wordno,
       bits;
  bitvector_t check,
              temp,
              carryxh,
              carryph,
              carrymh,
              Ph=0,
              Mh=0,
              Xv=0,
              Xh=0;

  res.a = -1;
  res.b = qlen;
  wordno = qlen/BITVECTOR_WORDSIZE;
  bits   = qlen & (BITVECTOR_WORDSIZE-1);
  wordno++;

  Pv  = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  Mv  = initbitvector(space, wordno*BITVECTOR_WORDSIZE);

  check = 0;
  bitvector_setbit(&check, bits, 1);
 
  setbitvector(Pv, wordno*BITVECTOR_WORDSIZE, 1);
  setbitvector(Mv, wordno*BITVECTOR_WORDSIZE, 0);


  for(i=0; i < slen; i++) {
   
    Eq = peq[enctab[(Uint)subject[i]]];
    carryxh = carryph = carrymh = 0;

    for(j=0; j < wordno; j++) {
     
      Xv = Eq[j] | Mv[j];
      temp = ((Eq[j] & Pv[j]) + Pv[j] + carryxh); 
      Xh = (temp ^ Pv[j]) | Eq[j];
      
      if (carryxh)
        carryxh = (temp <= (Eq[j] & Pv[j]) || temp <= Pv[j]);
      else
        carryxh = (temp < (Eq[j] & Pv[j]) || temp < Pv[j]);

      Ph = Mv[j] | ~(Xh | Pv[j]);
      Mh = Pv[j] & Xh;

      //check if last word
      if (j == wordno-1) {
        if (Ph & check) 
          score+=1;
        else if(Mh & check) 
          score-=1;
      }

      /*Ph = Ph << 1; with carry*/
      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp;

      Pv[j] = Mh | ~(Xv | Ph);
      Mv[j] = Ph & Xv;

    }

    if (score <= k && score <= res.b) {
        res.a = i;
        res.b = score;
    } 
  }
 
  FREEMEMORY(space, Pv);
  FREEMEMORY(space, Mv);       

  return res;
}



/*------------------------------ myersbitmatrix ------------------------------
 *    
 * @brief modified bitvector algorithm to return bitmatrix for backtracking
 * @author Steve Hoffmann 
 *   
 */
 
bitvector*
myersbitmatrix(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim) {

  bitvector 
  *Pv,
  *Mv,
  MvP,
  PvP,
  Eq;

  Uint score=qlen,
       i,
       j,
       wordno,
       bits;
  bitvector_t check,
              temp,
              carryxh,
              carryph,
              carrymh,
              Ph=0,
              Mh=0,
              Xv=0,
              Xh=0;

  res->a = -1;
  res->b = qlen;
  wordno = qlen/BITVECTOR_WORDSIZE;
  bits   = qlen & (BITVECTOR_WORDSIZE-1);
  wordno++;

  Pv = D; 
  Mv = &Pv[dim+1];
 
  memset(Pv[0], 255, wordno*(sizeof(bitvector_t)));
  memset(Mv[0], 0, wordno*(sizeof(bitvector_t)));

  check = 0;
  bitvector_setbit(&check, bits, 1);
 

  for(i=0; i < slen; i++) {
   
    Eq = peq[enctab[(Uint)subject[i]]];
    carryxh = carryph = carrymh = 0;

    MvP = Mv[i];
    PvP = Pv[i];

    for(j=0; j < wordno; j++) {
     
      Xv = Eq[j] | MvP[j];
      temp = ((Eq[j] & PvP[j]) + PvP[j] + carryxh); 
      Xh = (temp ^ PvP[j]) | Eq[j];
      
      if (carryxh)
        carryxh = (temp <= (Eq[j] & PvP[j]) || temp <= PvP[j]);
      else
        carryxh = (temp < (Eq[j] & PvP[j]) || temp < PvP[j]); 

      Ph = MvP[j] | ~(Xh | PvP[j]);
      Mh = PvP[j] & Xh;

      //check if last word
      if (j == wordno-1) {
        if (Ph & check) 
          score+=1;
        else if(Mh & check) 
          score-=1;
      }

      /*Ph = Ph << 1; with carry*/
      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp;

      Pv[i+1][j] = Mh | ~(Xv | Ph);
      Mv[i+1][j] = Ph & Xv;

    }

    if (score <= k && score <= res->b) {
        res->a = i;
        res->b = score;
    } 
  }

  return Pv;
}


/*---------------------------- bitvectorbacktrack ----------------------------
 *    
 * @brief backtracking in bitmatrix
 * @author Steve Hoffmann 
 *   
 */
 
Alignment*
bitvectorbacktrack(Alignment *al, bitvector *D, Uint dim, Uint k, Uint l) {
  int i=k-1, 
  j=l;
  bitvector *Pv = D;
  bitvector *Mv = &D[dim+1];
 
  while (i > 0 && j > 0) {
    if (bitvector_getbit(Pv[j], i)) {
      insertEop(al, Deletion);
      i--;
    } else {
      if (bitvector_getbit(Mv[j-1], i)) {   
        insertEop(al, Insertion);
      } else {
        insertEop(al, Replacement);
        i--;
      }
      j--;
    }
  }

  if (i==0 && j > 0) {
    insertEop(al, Replacement); 
  } else {
     while(i>=0) {
      insertEop(al, Deletion);
      i--;
    }
    /*no insertions in app. string matching at the end*/ 
  }

  /*adjust subject boundaries*/
  if(j>0) {
    al->voff = j-1;
    al->vlen -= (al->vlen > j) ?  j : al->vlen;
  }

  revMeops(al);
  return al;
}

/*------------------------------ wrapBitmatrix -------------------------------
 *    
 * @brief destruct bit matrix
 * @author Steve Hoffmann 
 *   
 */
 
void
wrapBitmatrix(void *space, bitvector *D, Uint m) {
  Uint i;
  for(i=0; i < m; i++) {
    FREEMEMORY(space, D[i]);
  }

  return;
}

