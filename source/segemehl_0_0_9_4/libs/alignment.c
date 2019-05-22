
/*
 *  alignment.c
 *  implementation to handle alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 02/03/2009 02:50:06 PM CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "mathematics.h"
#include "basic-types.h"
#include "alignment.h"
const char decodeEop[] = {'R','D','I'};

void 
initAlignment(Alignment *al, 
    char *u, Uint ulen, Uint uoff, 
    char *v, Uint vlen, Uint voff) {
  
  assert(uoff < ulen && voff < vlen);
  
  al->u = u;
  al->v = v;
  al->ulen = ulen;
  al->vlen = vlen;
  al->uoff = uoff;
  al->voff = voff;
  al->numofmeops = 0;
  al->meops = malloc(sizeof(Multieop)*(ulen+vlen));
  memset(al->meops, 0, sizeof(Multieop)*(ulen+vlen));
}

void 
wrapAlignment(Alignment *al) {
  free(al->meops);
  al->numofmeops = 0;
  al->u = NULL;
  al->v = NULL;
  al->vlen = 0;
  al->ulen = 0;
  al->uoff = 0;
  al->voff = 0;
}

void
countEops(Alignment *al, Uint *mat, Uint *mis, Uint *ins, Uint *del) {
  Uint i,j,k=0,l=0;
  
  *mat = 0;
  *mis = 0;
  *ins = 0;
  *del = 0;

  for(i=0; i < (*al).numofmeops; i++) {
    if((*al).meops[i].eop != Replacement) {
      if ((*al).meops[i].eop == Deletion) {
        *del += (*al).meops[i].steps;
        k += (*al).meops[i].steps;
      } else {
        *ins += (*al).meops[i].steps;
        l += (*al).meops[i].steps;
      }
    } else {
      for(j=0; j < (*al).meops[i].steps; j++) {
        if((*al).u[k+(*al).uoff] != (*al).v[l+(*al).voff]) 
          *mis += 1; 
        else 
          *mat += 1;
        
        k++; l++;
      }
    }
  }
  return;
}

Uint
getEdist(Alignment *al) {
  Uint i,j,k=0,l=0;
  Uint edist=0;
  for(i=0; i < (*al).numofmeops; i++) {
    if((*al).meops[i].eop != Replacement) {
      edist += (*al).meops[i].steps;
      if ((*al).meops[i].eop == Deletion) {
        k+= (*al).meops[i].steps;
      } else {
        l+= (*al).meops[i].steps;
      }
    } else {
      for(j=0; j < (*al).meops[i].steps; j++) {
        if((*al).u[k+(*al).uoff] != (*al).v[l+(*al).voff]) edist++;
        k++; l++;
      }
    }
  }
  return edist;
}

//dumps visual representation of alignments and should be shorter!
void 
showmultieoplist(Alignment *al) {

  Uint i;
  printf("[");
  for(i=0; i < (*al).numofmeops-1; i++) {
    printf("%c %d, ", decodeEop[(*al).meops[i].eop], (*al).meops[i].steps);
  }
  printf("%c %d]\n",decodeEop[(*al).meops[i].eop], (*al).meops[i].steps);    
}

//dumps visual representation of alignments and should be shorter!
char *
multieopstring(Alignment *al) {
  Uint i,j,q=0, p=0, cur=0, strsize, steps, msteps, ssteps;
  char *meopstr;
  char eopc=0;

  meopstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen)+1));

  for(i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    steps=0;
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      msteps=0;
      ssteps=0;
      for (j=0; j < al->meops[i].steps; j++) {
        if (al->u[j+p+al->uoff] != al->v[j+q+al->voff]) {
          if (j==0 || eopc == 'S') {
            ssteps++;
          } else {
            strsize = floor(log(msteps)/log(10))+3;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", msteps);
            cur+=strsize;
            msteps=0;
            ssteps=1;
          }
          eopc = 'S';
        } else {
          if (j==0 || eopc == 'M') {
            msteps++;
          } else {
            strsize = floor(log(ssteps)/log(10))+3;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", ssteps);
            cur+=strsize;
            msteps=1;
            ssteps=0;
          }
          eopc = 'M';
        }
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    } 
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      p+=steps;
    }
    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      eopc = 'I';  
      steps = al->meops[i].steps;
      q+=steps;
    }

    strsize = floor(log(steps)/log(10))+3;
    meopstr[cur] = eopc;
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }

  return meopstr;
}


//shows muliteoplist of all Alignments in *al
void 
showDynmultieoplist(Alignment* al, int size) {

  int i;
  for (i=0; i < size; i++) {
    showmultieoplist(&al[i]);
  }
}

//dumps visual representation of alignments and should be shorter!
void 
showAlign(Alignment* al, FILE *dev) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (utemp[j+r] != vtemp[j+r])
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      p+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      r+=j;
      q+=j;
    }
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s\n", len, &utemp[k*60]);
        fprintf(dev, "%.*s\n", len, &comp[k*60]);
        fprintf(dev, "%.*s\n", len, &vtemp[k*60]);
      }
      fprintf(dev, "\n");
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}

void
insertEop(Alignment *al, Eoptype eop) {

  //if previous multieops have been set up
  if (al->numofmeops > 0)  { 
    //inc steps if curr eop matches prev eops
    if (al->meops[(*al).numofmeops-1].eop == eop) {
      al->meops[(*al).numofmeops-1].steps++;
      //set up new multieop otherwise
    } else {
      al->numofmeops++;
      al->meops[(*al).numofmeops-1].eop =  eop;
      al->meops[(*al).numofmeops-1].steps = 1;
    }
    //set up first multieop
  } else {
    al->numofmeops = 1;
    al->meops[0].eop = eop;
    al->meops[0].steps = 1;
  }
}


void
revMeops(Alignment *al) {
  Uint start = 0;
  Uint end = al->numofmeops-1;
  Multieop *meops = al->meops; 

  if (al->numofmeops == 0) return;

  while (start<end) {
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;
    meops[end].eop ^= meops[start].eop;
    meops[end].steps ^= meops[start].steps;
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;

    start++;
    end--;
  } 
}






