
/*
 *  kdmatch.c
 *  routines 4 relaxed alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/27/2007 04:08:39 PM CET
 *
 *  SVN
 *  Revision of last commit: $Rev: 103 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-10 15:18:18 +0100 (Wed, 10 Dec 2008) $
 *
 *  Id: $Id: kdmatch.c 103 2008-12-10 14:18:18Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/kdmatch.c $
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "memory.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include "manout.h"
#include "biofiles.h"
#include "vtprogressbar.h"
#include "karlin.h"
#include "sort.h"
#include "basic-types.h"
#include "bitvectoralg.h"
#include "bitVector.h"
#include "kdmatch.h"
#include "bitArray.h"
#include "segemehl.h"
#include "container.h"
#include "kdchain.h"
#include "debug.h"
#include "info.h"
#include "kdseed.h"
#include "alignment.h"
#include "sw.h"
#include "seqclip.h"
#include <pthread.h>
#include "iupac.h"




/*-------------------------- bl_compareSliceEvents ---------------------------
 *    
 * @brief compare function to sort splice events in read order
 * @author Steve Hoffmann 
 *   
 */

int
bl_compareSplits (const void *a, const void *b)
{
  split_t *eventa, *eventb;

  eventa = (split_t*) a; 
  eventb = (split_t*) b; 
	
  if(eventa->i > eventb->i) return 1;
  if(eventa->i < eventb->i) return  -1;

  return 0;
}


/*----------------------- bl_compareSpliceEventMapElem -----------------------
 *    
 * @brief compare the elems of the Splice Event Map Elems for remapping
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_compareSpliceEventMapElem(const void *a, const void *b )
{
    spliceeventmapelem_t *elema, *elemb;
    Uint idxa, idxb, posa, posb;

    elema = (spliceeventmapelem_t*) a;
    elemb = (spliceeventmapelem_t*) b;

    idxa = elema->ptr->subidx[elema->site];
    idxb = elemb->ptr->subidx[elemb->site];

    if(idxa > idxb) return 1;
    if(idxa < idxb) return -1;

    if(elema->type == 0) { 
      posa = elema->ptr->start[elema->site];
    } else { 
      posa = elema->ptr->end[elema->site];
    }

    if(elemb->type == 0) { 
      posb = elemb->ptr->start[elemb->site];
    } else {
      posb = elemb->ptr->end[elemb->site];
    }

    if(posa > posb) return 1;
    if(posb < posb) return -1;

 
	return 0;
}


/*------------------------- bl_destructSpliceEvents --------------------------
 *    
 * @brief destruct the list of splice events
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_destructSpliceEvents (spliceevents_t *events)
{
  Uint i;

  for(i=0; i < events->noofevents; i++) {
    if(events->event[i].subidx) { 
      FREEMEMORY(space, events->event[i].subidx);
      events->event[i].subidx = NULL;
    }
    if(events->event[i].start) {
       FREEMEMORY(space, events->event[i].start);
       events->event[i].start = NULL;
    }
    if(events->event[i].end) {
      FREEMEMORY(space, events->event[i].end);
      events->event[i].end = NULL;
    }
    if(events->event[i].strand) {
      FREEMEMORY(space, events->event[i].strand);
      events->event[i].strand = NULL;
    }
    
  }

  FREEMEMORY(space, events->event);
  return ;
}


/*---------------------- bl_reconstructSpliceTranscript ----------------------
 *    
 * @brief reconstruct a spliced transcript from splice event
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_reconstructSpliceTranscript (void *space, spliceevent_t *event, MultiCharSeq *seq, 
    char **seqs, Uint off)
{
  
  Uint i, refstart, refend, reflen, chrstart, chrend, curlen = 0,
       starti, startu=0, startp=0, endj = 0, endu=0, endq=0;
  Uint fiveprimestart, threeprimestart, fiveprimelen, threeprimelen;
  char *temp=NULL, *refseq, *spliceseq=NULL, *fiveprime, *threeprime;


  starti = strlen(seqs[0]);
  for(i=0; i < event->noofsites; i++) {
      
    getMultiCharSeqIdxBounds(seq, event->subidx[i], &chrstart, &chrend);
      assert(chrstart <= event->start[i]);
      assert(chrend >= event->end[i]);

      refstart = event->start[i];
      refend = event->end[i];

      fiveprime = NULL;
      fiveprimelen = 0;

    if(i == 0) {

      if(event->strand[i] ==1) {
        fiveprimestart = (chrend - refend > off) ? refend + off : chrend - 1;
        if(fiveprimestart > refend+1) { 
          fiveprimelen = fiveprimestart-refend;
          fiveprime = charIUPACcomplement(space, &seq->sequences[refend+1], fiveprimelen);
        }
      } else {
        fiveprimestart = (refstart - chrstart > off) ? refstart - off : chrstart + 1;
        if(fiveprimestart+1 < refstart) {
          fiveprime = ALLOCMEMORY(space, NULL, char, refstart-fiveprimestart+1);
          fiveprimelen = refstart-fiveprimestart;
          memmove(fiveprime, &seq->sequences[fiveprimestart], fiveprimelen);
          fiveprime[fiveprimelen] = 0;
        }
      }

      fprintf(stdout, "5':%d (%d)\n", fiveprimestart, fiveprimelen);
      spliceseq = ALLOCMEMORY(space, spliceseq, char, curlen+fiveprimelen+1);
      memmove(&spliceseq[curlen], fiveprime, fiveprimelen);
      curlen += fiveprimelen;
    }

    reflen = (refend - refstart) +1;
    refseq = &seq->sequences[refstart];

    if(event->i[i] < starti) { 
      startu = event->strand[i];
      startp = event->start[i];
      starti = event->i[i];
    }

    if(event->j[i] > endj) {
      endu = event->strand[i];
      endq = event->end[i];
      endj = event->j[i];
    }

    if(event->strand[i] ==1) {
      temp = charIUPACcomplement(space, refseq, reflen);
      fprintf(stdout, "- :");
    } else {
      fprintf(stdout, "+ :");
      temp = ALLOCMEMORY(space, NULL, char, reflen+1);
      memmove(temp, refseq, reflen);
      temp[reflen] = 0;
    }

    fprintf(stdout, "[%d,%d]->[%d,%d]\n", event->i[i], event->j[i], event->start[i], event->end[i]);

    spliceseq = ALLOCMEMORY(space, spliceseq, char, curlen+reflen+1);
    memmove(&spliceseq[curlen], temp, reflen);
    curlen += reflen;

    if(i == event->noofsites-1) {
      threeprime = NULL;
      threeprimelen = 0;

      if(event->strand[i]==1) {
        threeprimestart = (refstart - chrstart > off) ? refstart - off : chrstart + 1;
        if(threeprimestart+1 < refstart) {
          threeprime = ALLOCMEMORY(space, NULL, char, refstart-threeprimestart+1);
          threeprimelen = refstart-threeprimestart;
          memmove(threeprime, &seq->sequences[threeprimestart], threeprimelen);
          threeprime[threeprimelen] = 0;
        }
      } else {
        threeprimestart = (chrend - refend > off) ? refend + off : chrend - 1;
        if(threeprimestart > refend+1) { 
          threeprimelen = threeprimestart-refend;
          threeprime = charIUPACcomplement(space, &seq->sequences[refend+1], threeprimelen);
        }
      }

      fprintf(stdout, "3':%d (%d)\n", threeprimestart, threeprimelen);
      spliceseq = ALLOCMEMORY(space, spliceseq, char, curlen+threeprimelen+1);
      memmove(&spliceseq[curlen], threeprime, threeprimelen);
      curlen += threeprimelen;
    }

    FREEMEMORY(space, temp);
  }

  fprintf(stdout, "\n");
  
  if(spliceseq) { 
    spliceseq[curlen] = 0;

    fprintf(stdout, "0:%s\n", seqs[0]);
    fprintf(stdout, "1:%s\n", seqs[1]);
    fprintf(stdout, "S:%s\n", spliceseq);
    fprintf(stdout, "(%d,%d) -> (%d,%d)\n", starti, startu, endj, endu);
    FREEMEMORY(space, spliceseq);
  }

  return ;
}

/*--------------------------- bl_storeSpliceEvent ----------------------------
 *    
 * @brief store the splice event for remapping
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_storeSpliceEvent (void *space, MultiCharSeq *seq, 
    gmatchlist_t *list, spliceevents_t *events, Uint readid, 
    Uint off, char **seqs, char *name)
{
  Uint u, i, noofsplits=0;
  gmatch_t *match =NULL;
  spliceevent_t *event;
  split_t *splits = NULL;
  
   

  for(u=0; u < 2; u++) {
    for(i=0; i<list->n[u]; i++) {

      match = &list->matches[u][i];
      splits = ALLOCMEMORY(space, splits, split_t, noofsplits+1);
      splits[noofsplits].subidx = match->subject;
      splits[noofsplits].start=  match->p;
      splits[noofsplits].end = match->q;
      splits[noofsplits].strand = u;
      splits[noofsplits].i = match->i;
      splits[noofsplits].j = match->j;
      noofsplits++;
    }
  }

  qsort(splits, noofsplits, sizeof(split_t), bl_compareSplits);

  event = ALLOCMEMORY(space, NULL, spliceevent_t, 1);
  event->firstreadid = readid;
  event->noofsites = 0;
  event->subidx = NULL;
  event->start = NULL;
  event->end = NULL;
  event->strand = NULL;
  event->i = NULL;
  event->j = NULL;
  fprintf(stdout, "name: %s\n", name);

  for(i=0; i< noofsplits; i++) {
    event->subidx = ALLOCMEMORY(space, event->subidx, Uint, event->noofsites+1);
    event->start = ALLOCMEMORY(space, event->start, Uint, event->noofsites+1);
    event->end = ALLOCMEMORY(space, event->end, Uint, event->noofsites+1);
    event->strand = ALLOCMEMORY(space, event->strand, char, event->noofsites+1);
    event->i = ALLOCMEMORY(space, event->i, uint16_t, event->noofsites+1);
    event->j = ALLOCMEMORY(space, event->j, uint16_t, event->noofsites+1);
    event->subidx[event->noofsites] = splits[i].subidx;
    event->start[event->noofsites] = splits[i].start;
    event->end[event->noofsites] = splits[i].end;
    event->strand[event->noofsites] = splits[i].strand;
    event->i[event->noofsites] = splits[i].i;
    event->j[event->noofsites] = splits[i].j;
    event->noofsites++;
  }

  bl_reconstructSpliceTranscript (space, event, seq, seqs, off);
  FREEMEMORY(space, splits);
  //events->noofevents++;
  return ;
}

/*----------------------- bl_kdConstructSpliceEventMap -----------------------
 *    
 * @brief splice map construction
 * @author Steve Hoffmann 
 *   
 */

spliceeventmap_t *
bl_kdConstructSpliceEventMap (spliceevents_t *events)
{
  Uint inc = 1000;
  Uint i, j, size = 1000, n = 0;
  spliceeventmap_t *eventmap;
  spliceeventmapelem_t *map;
  spliceevent_t *event;

  eventmap = ALLOCMEMORY(space, NULL, spliceeventmap_t, 1);
  map = ALLOCMEMORY(space, NULL, spliceevents_t*, size);

  for(i=0; i < events->noofevents; i++) {
    event = &events->event[i];
    for(j=0; j < event->noofsites; j++) {
               
      if(n+2 >= size) {
        size += inc;
        map = ALLOCMEMORY(space, map, spliceeventmapelem_t, size);
      }

      if(j > 0) { 
        map[n].site = j;
        map[n].ptr = event;
        map[n].type = 0;
        n++;
      }
      if(j < event->noofsites-1) {
        map[n].site = j;
        map[n].ptr = event;
        map[n].type = 1;
        n++;
      }
    }
  }

  qsort(map, n, sizeof(spliceeventmapelem_t), bl_compareSpliceEventMapElem);

  eventmap->map = map;
  eventmap->size = n;

  return eventmap;
}




/*-------------------------- bl_kdRemapGetBestSeed ---------------------------
 *    
 * @brief get the best seed from fasta description
 * @author Steve Hoffmann 
 *   
 */
 
bestseed_t*
bl_kdRemapGetBestSeed (void *space, char *desc, Uint desclen)
{
  char *copy, *cur, strand;
  Uint clen, idx, pos, mat, start;
  bestseed_t* seed;

  cur = desc + (desclen-1);
  while(*cur != ';') {
    cur--;
  }
  
  clen = strlen(cur);
  copy = ALLOCMEMORY(space, NULL, char, clen+1);
  memmove(copy, cur, clen);
  copy[clen] = 0;


  /* read readstart */
  cur = strtok(copy, ":");
  if(cur == NULL) return NULL;
  start = atoi(cur);
  if(cur[0] != '0' && start == 0) return NULL;
  

  /* read seedlen */
  cur = strtok(NULL, ":");
  if(cur == NULL) return NULL;
  mat = atoi(cur); 
  if(cur[0] != '0' && mat == 0) return NULL;
 

  /* read idx */
  cur = strtok(NULL,":");
  if(cur == NULL) return NULL;
  idx = atoi(cur);
  if(cur[0] != '0' && idx == 0) return NULL;

  /* read refpos */
  cur = strtok(NULL, ":");
  if(cur == NULL) return NULL;
  pos = atoi(cur);
  if(cur[0] != '0' && pos == 0) return NULL;
  
  /* read strand */
  cur = strtok(NULL, ":");
  if(cur == NULL) return NULL;
  strand = cur[0];
  if(cur[0] != '+' && cur[0] != '-') return NULL;

  FREEMEMORY(space, copy);

  seed = ALLOCMEMORY(space, NULL, bestseed_t, 1);
  seed->readstart = start;
  seed->mat = mat;
  seed->refidx = idx;
  seed->refpos = pos;
  seed->refstrand = strand;


  return seed;
}


/*-------------------------------- bl_kdRemap --------------------------------
 *    
 * @brief remapping the reads to splice sites
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_kdRemap (void *space, spliceeventmap_t  *map, fasta_t *reads)
{
  Uint k, desclen;
  bestseed_t *seed;
  char *desc;
	 
  for (k=0; k < reads->noofseqs; k++) {
    desc = bl_fastaGetDescription(reads, k);
    desclen = bl_fastaGetDescriptionLength(reads, k);
    seed = bl_kdRemapGetBestSeed(space, desc, desclen);

    if(seed) FREEMEMORY(space, seed);
  }

  return ;
}



/*---------------------------- se_kdFindBestMate -----------------------------
 *    
 * @brief find the 'best' pair from two lists of hits for query and mate
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdFindBestMatePair (void *space, gmatchlist_t *querylist, 
    gmatchlist_t *matelist) {

  Uint u, v, i, j, ucnt=0, vcnt=0, uprime=0, vprime=0;
  unsigned char found = 0, downstream = 0;
  Alignment *al;
  PairUint p,q;
  gmatchlist_t *list;

  list = NULL;
  memset(&p, 0, sizeof(PairUint));
  memset(&q, 0, sizeof(PairUint));

  for(u=0; u < 2; u++) {
    ucnt += querylist->n[u];
    for(i=0; i < querylist->n[u]; i++) {
      uprime = u;
      for(v=0; v < 2; v++) {
        vcnt += matelist->n[v];
        for(j=0; j < matelist->n[v]; j++) {
          vprime = v;
          if(querylist->matches[u][i].subject == 
              matelist->matches[v][j].subject) {
            if(!found || abs((LLint)querylist->matches[p.a][p.b].p - 
                  matelist->matches[q.a][q.b].p) >
                abs((LLint)querylist->matches[u][i].p -
                  matelist->matches[v][j].p)) 
            {
              found = 1;
              p.a = u;
              p.b = i;
              q.a = v;
              q.b = j;
            }
          }
        }
      }
    }
  }


  if(!found && ucnt == 1 && vcnt == 1) { 
    p.a = uprime;
    p.b = 0;
    q.a = vprime;
    q.b = 0;
    found=1;
  }


  if(found) {

    al = ALLOCMEMORY(space, NULL, Alignment, 1);
    copyAlignment(al, querylist->matches[p.a][p.b].al);
   
    list = bl_gmatchlistInit(space, querylist->matches[p.a][p.b].edist);  
    list = se_kdMatchListAdd(list,  
        querylist->matches[p.a][p.b].subject,
        querylist->matches[p.a][p.b].p,
        querylist->matches[p.a][p.b].q,
        querylist->matches[p.a][p.b].edist,
        querylist->matches[p.a][p.b].scr,
        querylist->matches[p.a][p.b].i,
        querylist->matches[p.a][p.b].j-1,
        querylist->matches[p.a][p.b].evalue, al, p.a, -1, -1, 0, -1, -1, 0, 0);
        
    if(querylist->matches[p.a][p.b].p >= matelist->matches[q.a][q.b].p) 
      downstream = 1; 
    else 
      downstream = 0; 

    al = ALLOCMEMORY(space, NULL, Alignment, 1);
    copyAlignment(al, matelist->matches[q.a][q.b].al);

    se_kdSetMate(space, &list->matches[p.a][0], 
//        list->matches[p.a][0].subject, 
        matelist->matches[q.a][q.b].subject,
        matelist->matches[q.a][q.b].p, 
        matelist->matches[q.a][q.b].q, 
        matelist->matches[q.a][q.b].edist, 
        al, downstream, (p.a != q.a));

    if(list->mateminedist > matelist->matches[q.a][q.b].edist) {
      list->mateminedist = matelist->matches[q.a][q.b].edist;
    }

    querylist->matches[p.a][p.b].skip = 1;
    matelist->matches[q.a][q.b].skip = 1; 
  }

  return list;
}


/*------------------------------ se_kdAlignMate ------------------------------
 *    
 * @brief find the mate once a sequence was located
 * @author Steve Hoffmann 
 *   
 */


  Uint
se_kdAlignMate(void *space, MultiCharSeq *seq, char **seqs, Uint len, 
    gmatchlist_t *list, Uint maxedist,Uint* enctab, bitvector *D, Uint maxlen) 
{

  PairSint mb;
  Alignment *al;
  bitvector *peq[2];
  char *refseq, *upstreamrefseq;
  Uint u, i, k, p;
  Uint idx, refstart, reflen, upstreamreflen, 
       upstreamrefstart, chrstart, chrend;

  peq[0] = getpeq(space, seqs[0], len, seq->map, 
      seq->mapsize, enctab);
  peq[1] = getpeq(space, seqs[1], len, seq->map, 
      seq->mapsize, enctab);
 

  for(u=0; u < 2; u++) {
    for (i=0; i < list->n[u]; i++) {

      idx = list->matches[u][i].subject;
      getMultiCharSeqIdxBounds(seq, idx, &chrstart, &chrend);

      p = list->matches[u][i].p;
      refstart = p;
      reflen = (chrend > (Lint)refstart + maxlen)? maxlen :(chrend-refstart);
      refseq = &seq->sequences[refstart];

      upstreamreflen =((Lint)p-maxlen > chrstart)? maxlen :(Lint)p -chrstart;
      upstreamrefstart = p - upstreamreflen;
      upstreamrefseq = &seq->sequences[upstreamrefstart];

      for(k=0; k < 2; k++) {

        myersbitmatrix(NULL, seqs[k], len, refseq, reflen, 
            seq->map, seq->mapsize, enctab, 
            len-maxedist, peq[k], &mb, D, reflen);
          

        if (mb.a != -1 && mb.b <= maxedist && mb.a < reflen) {  
          al = ALLOCMEMORY(space, NULL, Alignment, 1);

          initAlignment(al, seqs[k], len, 0, refseq, reflen, 0);
          bitvectorbacktrack(al, D, reflen, len, mb.a);

          se_kdSetMate(space, &list->matches[u][i], idx, 
              refstart+al->voff, refstart+mb.a-1, 
              mb.b, al, 1, (u != k));

          if(list->mateminedist > mb.b) {
            list->mateminedist = mb.b;
          }
          

          list->matches[u][i].noofmatematches++;
        }

        myersbitmatrix(NULL, seqs[k], len, upstreamrefseq, upstreamreflen, 
            seq->map, seq->mapsize, enctab, 
            len-maxedist, peq[k], &mb, D, upstreamreflen);

        
        if (mb.a != -1 && mb.b <= maxedist && mb.a < upstreamreflen) {  
          al = ALLOCMEMORY(space, NULL, Alignment, 1);

          initAlignment(al, seqs[k], len, 0, upstreamrefseq, 
              upstreamreflen, 0);
          bitvectorbacktrack(al, D, upstreamreflen, len, mb.a);

          se_kdSetMate(space, &list->matches[u][i], idx, 
              upstreamrefstart+al->voff, upstreamrefstart+mb.a-1, 
              mb.b, al, 0, (u != k));

          if(list->mateminedist > mb.b) {
            list->mateminedist = mb.b;
          }
  

          list->matches[u][i].noofmatematches++;
        }
      }
    }
  }

  for(u=0; u < 2; u++) {
    for(i=0; i < seq->mapsize; i++) {
      FREEMEMORY(space, peq[u][i]);
    }  
    FREEMEMORY(space, peq[u]);
  }

  return 0;
}



/*--------------------------- bl_kdUpdateBestSeed ----------------------------
 *    
 * @brief update the best seed record
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_kdUpdateBestSeed (bestseed_t *best, MultiCharSeq *seq, Suffixarray *s, 
    Uint readstart, Uint mat, Uint l, Uint refstrand)
{
  
  Uint pos, subidx, substart, subend;    
  
  pos = s->suftab[l];
  subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, subidx, &substart, &subend);
  pos -= substart;

  if(best->mat < mat) {
    best->readstart = readstart;
    best->mat= mat;
    best->refidx = subidx;
    best->refpos = pos;
    best->refstrand = refstrand;
  }
	
  return ;
}

/*--------------------------- se_kdMatchStemAlign ----------------------------
 *    
 * @brief align the seeds in the matchstem using bv algorithmics
 * @author Steve Hoffmann 
 *   
 */


gmatchlist_t*
se_kdMatchStemAlign(void *space, Suffixarray *s, MultiCharSeq *seq, 
    matchstem_t **stems, char **seqs, Uint len, karlin_t *stats, 
    segemehl_t *nfo, Uint *enctab, bitvector* D, bestseed_t *best) {

  Uint u, k, j, l, r, q, i, pos, mat, mis, ins, del; 
  int maxedist, bestedist, scr, skipmargin=0;
  Uint *check=NULL;
  Uint checklen=0;
  double E;
  bitvector *peq[2];
  PairSint mb;
  MultiCharSeqAlignment mcsa;
  gmatchlist_t *list;

  maxedist = bestedist = len-ceil((nfo->accuracy*len)/100);
  skipmargin = 40*((double)maxedist/100.);
  list = bl_gmatchlistInit(space, maxedist);

  peq[0] = getpeq(NULL, seqs[0], len, seq->map, 
      seq->mapsize, enctab);
  peq[1] = getpeq(NULL, seqs[1], len, seq->map, 
      seq->mapsize, enctab);

  for(u = 0; u < 2; u++) {
    for(i = 0; i < len; i++) {
      for(q = 0; q < stems[u][i].noofbranches; q++) {

        l = stems[u][i].branches[q].l; 
        r = stems[u][i].branches[q].r;
        mat = stems[u][i].branches[q].mat;
        mis = stems[u][i].branches[q].mis;
        ins = stems[u][i].branches[q].ins;
        del = stems[u][i].branches[q].del;

        E = kd_getBranchEvalue(stems[u], i, q, len, s->numofsuffixes, stats);

        if (l > r || E > nfo->maxevalue || (r-l) > nfo->M) 
          continue;
          
        bl_kdUpdateBestSeed(best, seq, s, i, mat, l, u);

        for(j = l; j <= r; j++) {
          pos = s->suftab[j];

          if(mat != len || mis+ins+del != 0) {
            initMultiCharSeqAlignment(space, &mcsa, seq, pos, 
                i+maxedist, len+2*(maxedist+1), u, NULL, seqs[u], len);
          } else {
            initMultiCharSeqAlignment(space, &mcsa, seq, pos, 
                0, len, u, NULL, seqs[u], len);
          }

          /*skip or update identical matches*/
          for(k = 0; k < checklen; k++) 
            if (check[k] >= mcsa.refstart-skipmargin && 
                check[k] <= mcsa.refstart+skipmargin)     
              break;	

          if (k < checklen) {
            wrapMultiCharSeqAlignment(space, &mcsa);
            continue;
          }

          check = ALLOCMEMORY(space, check, Uint, checklen+1);
          check[checklen++]= mcsa.refstart;

          if (mat == len && mis+ins+del == 0) {
            scr = kd_getBranchScore(stems[u], i, q);

            for(k=0; k < len; k++) {
              insertEop(mcsa.al, Replacement);
            }

            mb.b = getEdist(mcsa.al);

            if (mb.b <= maxedist && mb.b <= bestedist){	      
              
              list = se_kdMatchListAdd(list, mcsa.subidx, 
                  pos, pos+len-1, mb.b, scr, 0, len-1, E, mcsa.al, 
                  u, -1, -1, 0, -1, -1, 0, 0);

              if(nfo->bestonly) bestedist = list->minedist;
            
            } else {
              wrapMultiCharSeqAlignment(space, &mcsa);
            }
          
          } else {

            myersbitmatrix(NULL, seqs[u], len, mcsa.refseq, mcsa.reflen, 
                seq->map, seq->mapsize, enctab, len-bestedist, peq[u], 
                &mb, D, mcsa.reflen);

            if (mb.a != -1 && mb.b <= maxedist && 
                mb.b <= bestedist && mb.a < mcsa.reflen) {  
              bitvectorbacktrack(mcsa.al, D, mcsa.reflen, len, mb.a);

              /*skip or update identical matches*/
              for(k = 0; k < list->n[u]; k++) 
                if (list->matches[u][k].p == mcsa.refstart+mcsa.al->voff)
                  break;	

              if (k < list->n[u]) { 
                if (list->matches[u][k].edist <= mb.b){
                  wrapMultiCharSeqAlignment(space, &mcsa);
                } else {
                  scr = kd_getBranchScore(stems[u], i, q);

                  list = se_kdMatchListSet(space, list, mcsa.subidx, 
                      mcsa.refstart+mcsa.al->voff, 
                      mcsa.refstart+mb.a-1,
                      mb.b, scr, 0, len-1, E, mcsa.al, u, k);
                }
                continue;
              }

              scr = kd_getBranchScore(stems[u], i, q); 

              list=se_kdMatchListAdd(list, mcsa.subidx, 
                  mcsa.refstart+mcsa.al->voff, 
                  mcsa.refstart+mb.a-1, mb.b, scr, 0, len-1, E, mcsa.al, 
                  u, -1, -1, 0, -1, -1, 0, 0);

              if(nfo->bestonly) bestedist = list->minedist;
            
            } else { 
              wrapMultiCharSeqAlignment(space, &mcsa);
            }
          }
        }
      }
    }

    for(j=0; j < seq->mapsize; j++) {
      FREEMEMORY(space, peq[u][j]);
    }

    FREEMEMORY(space, peq[u]);
    if(check) {
      FREEMEMORY(space, check);
      check = NULL;
      checklen = 0;
    }
  }

  return list;
}


/*--------------------------- se_kdAlignSplitChain ---------------------------
 *    
 * @brief align a chain of fragments using local multi spliced alignment
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdAlignSplitChain (void *space, branchChain_t *chains, Uint noofchains,
    Suffixarray *arr, MultiCharSeq *seq, char *querydesc, matchstem_t **stems,
    char **seqs, Uint qrylen, int *scores, int indel, spliceevents_t *events, segemehl_t *nfo) {

  Uint i, j, start, k, edist=0, prev = -1, fragno = 0,
    ulen=0, vlen=0, vstart=0, ustart=0, nextustart=0, prevustart=0, uend=0, maxedist, ustartj=0, 
    *strands, *reflens, totalcover = 0, totaledist = 0,
    q = 0, p =0, **bd, l, ll, r, rr, beststart=0, lastsubidx=0, d1idx=0, 
    d2idx=0, d3idx=0, sub_start, sub_end, previdx, prevpos, nextidx, nextpos;
  
  double mindist = DBL_MAX, d1=0, d2=0, di, dj;
  unsigned char rmalign=0, laststrand=0, trans=0, purge=0;
  char **refseqs, nextstrand='+', prevstrand='+';
  unsigned char bestchain = 0;
  int *M, **lmr, **lmv, **lmc, totalscore=0, score=0;
  branchChain_t *cur;
  MultiCharSeqAlignment *a, *b;
  Alignment **aligns, *alcopy; 
  gmatchlist_t *list=NULL;

  maxedist = qrylen - ceil((nfo->accuracy * qrylen)/100);
  list = bl_gmatchlistInit(space, maxedist);
 
  if(noofchains == 0) return list;

  qsort(chains, noofchains, sizeof(branchChain_t), cmp_chainscores);
  
  /*DEBUG SHOW CHAINS*/
  //showChains(chains, noofchains, arr, stdout, seqs[1], qrylen);

  cur = &chains[0];
  bestchain = 0;
  if (noofchains > 1) { 
    if (((double)chains[0].score)*0.8 > chains[1].score) {
      bestchain = 1;
    } else {
      return list;
    }
  }

  if(cur->nooffragments <= 1) return list; 

  l = cur->f[0]->branch->l;
  r = cur->f[0]->branch->r;

  bd = ALLOCMEMORY(space, NULL, Uint*, r-l+1);

  /* minimize distance of fragment hits within chain
   * for all possible start loci in [l,r] select
   * a chain of closest loci*/

  for(i=0, p=l; p <= r; p++, i++) {
    bd[i] = calloc(cur->nooffragments, sizeof(Uint));
    bd[i][0] = p;

    for(di=0, j=1; j < cur->nooffragments; j++) {
      ll = cur->f[j]->branch->l;
      rr = cur->f[j]->branch->r;

      for(dj=0, q=ll; q <= rr; q++) {
        d1idx = getMultiCharSeqIndex(seq, 
            &seq->sequences[arr->suftab[q]]);
        d2idx = getMultiCharSeqIndex(seq, 
            &seq->sequences[arr->suftab[bd[i][j-1]]]);

        if(d1idx != d2idx) {
          d1 = UINT_MAX;
        } else {
          d1 = llabs((LLint) arr->suftab[bd[i][j-1]] - arr->suftab[q]);
        }

        if(bd[i][j]) {
          d3idx = getMultiCharSeqIndex(seq, 
              &seq->sequences[arr->suftab[bd[i][j]]]);
          d2 = llabs((LLint) arr->suftab[bd[i][j-1]] - arr->suftab[bd[i][j]]);
        }

        if(d3idx != d2idx) {
          d2 = UINT_MAX;
        }
        
        if(!bd[i][j] || d1 < d2) {
          bd[i][j] = q;
          dj = d1;
        }       
      }
      di += dj;
    }

    if(di < mindist) {
      beststart = i;
      mindist = di;
    }
  }

  a = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, cur->nooffragments);
  b = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, cur->nooffragments);

  reflens = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  strands = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  refseqs = ALLOCMEMORY(space, NULL, char*, cur->nooffragments);
  aligns =  ALLOCMEMORY(space, NULL, Alignment*, cur->nooffragments);
  
  for(i=0; i < cur->nooffragments; i++) {
    
    start = arr->suftab[bd[beststart][i]]-20;
    initMultiCharSeqAlignment(space, &a[i], seq, start,
        cur->f[i]->start+maxedist, qrylen+2*(maxedist+1+20), cur->f[i]->strand,
        querydesc, seqs[cur->f[i]->strand], qrylen);
    
    aligns[i]  = a[i].al;
    refseqs[i] = a[i].refseq;
    reflens[i] = a[i].reflen; 
    strands[i] = cur->f[i]->strand;

//    fprintf(stdout, "%s\n fragment:%d, start:%d, strand:%d, curstart:%d, maxedist:%d\n", 
//        querydesc, i, start, strands[i],  cur->f[i]->start, maxedist);
//    fprintf(stdout, "\n");
    
  }

  M = localmultisplicedmatrix(space, seqs[0], seqs[1], qrylen,
      refseqs, reflens, strands, cur->nooffragments, indel,
      constscr, scores, &lmv, &lmr, &lmc);

  if(M == NULL) {
    fprintf(stderr, "empty matrix returned for seqs: '%s'/'%s' (%d)\n", 
        seqs[0], seqs[1], qrylen);

    for(i=0; i < cur->nooffragments; i++) {

      getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);
      fprintf(stderr, "fragment %d: %d in %d[%d,%d] '", 
          i, arr->suftab[bd[beststart][i]], a[i].subidx, sub_start, sub_end);
      for(j=0; j< qrylen; j++) fprintf(stderr, "%c", refseqs[i][j]);
      fprintf(stderr, "'(%d) strand:%d\n", reflens[i], strands[i]);
    }
    return list;
  }

  localmultisplicedtraceback(space, M, seqs[0], seqs[1], qrylen, 
      refseqs, reflens, strands, cur->nooffragments, indel,
      constscr, scores, aligns, lmv, lmr, lmc);

  purge = 0;
  for(k=0, i=0; i < cur->nooffragments; i++) {
  
    ustart = a[i].al->uoff;
    vstart = a[i].al->voff;
    ulen = getUalignlen(a[i].al);
    vlen = getValignlen(a[i].al);
    score = getSWScore(a[i].al, scores, indel);
    edist = getEdist(a[i].al);

  // fprintf(stdout, "%s\n ustart:%d, vstart:%d, strand:%d, score:%d, edist:%d\n", 
  //      querydesc, ustart, vstart, a[i].strand,  score, edist);

    if(edist > (ulen - ceil((nfo->accuracy * ulen)/100))) {  
  //    fprintf(stdout, "purging!\n");
      purge = 1;
    }
    
    if (ulen >= nfo->minfragmentalignlen && 
        vlen >= nfo->minfragmentalignlen &&
        score >= nfo->minfragmentalignscore) { 

      totalcover += ulen;
      totalscore += score;
      totaledist += edist;
      memmove(&b[k], &a[i], sizeof(MultiCharSeqAlignment));
      k++;
      
      alcopy = ALLOCMEMORY(space, NULL, Alignment, 1);
      copyAlignment(alcopy, a[i].al);     
      if (a[i].strand == 1) {
        uend = qrylen - ustart - 1;
        ustart = uend - ulen + 1;        
      } else {
        uend = ustart + ulen - 1;
      }

      previdx = -1;
      prevpos = -1;
      nextidx = -1;
      nextpos = -1;
      prevstrand = -1;
      nextstrand = -1;
      prevustart =0;
      nextustart = 0;


      for(j=0; j < cur->nooffragments; j++) {

        if(a[j].strand == 1) {
          ustartj = qrylen - a[j].al->uoff - getUalignlen(a[j].al);
        } else {
          ustartj = a[j].al->uoff;
        }


        if (ustartj < ustart &&  (!prevustart || ustartj >= prevustart) &&
            getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getSWScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 
          
          previdx = a[j].subidx;

          if(a[j].strand == 0) { 
            prevpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
            prevstrand = '+';
            if(a[j].strand == 1) { 
              prevstrand = '-';
            }
          } else {  
            prevpos = a[j].refstart + a[j].al->voff;
            prevstrand = '-';
            if(a[j].strand == 0) { 
              prevstrand = '+';
            }
          }
          prevustart = ustartj;
        }
      }

    
      for(j=0; j < cur->nooffragments; j++) { 

        if(a[j].strand == 1) {
          ustartj = qrylen - a[j].al->uoff - getUalignlen(a[j].al);
        } else {
          ustartj = a[j].al->uoff;
        }

        if (ustartj > ustart && (!nextustart || ustartj <= nextustart) &&
            getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getSWScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 
          nextidx = a[j].subidx;
          if(a[j].strand == 0) { 
            nextpos = a[j].refstart + a[j].al->voff;
            nextstrand = '+';
          } else { 
            nextpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
            nextstrand = '-';
          }
          nextustart = ustartj;
        }
      }

      list = se_kdMatchListAdd(list, a[i].subidx, 
          a[i].refstart + vstart, 
          a[i].refstart + vstart + vlen - 1, 
          edist, score, ustart, //ustart + ulen - 1, 
          uend, .0, alcopy, a[i].strand, 
          previdx, prevpos, prevstrand, 
          nextidx, nextpos, nextstrand, fragno);

      prev = i;
      fragno++;

      if (k > 1 && (laststrand != a[i].strand || lastsubidx != a[i].subidx)) {
        trans = 1;
      }

      laststrand = a[i].strand;
      lastsubidx = a[i].subidx;
    } 
  }

  totalcover *= 100;
  totalcover /= qrylen;

//  fprintf(stdout, "totalcover %d, totalscore %d, k %d, noofchains %d\n", 
//      totalcover, totalscore, k, noofchains);
  
  if(totalscore >= nfo->minsplicedalignscore && 
     totalcover >= nfo->minsplicedaligncover && k > 1 &&
     (!trans || (noofchains == 1 || bestchain)) && !purge) {
    /*restrictive policy for reporting trans splicing events*/
//    if(!trans || noofchains == 1) {
//      reportSplicedMatch(space, querydesc, b, k, 
//        totalcover, totaledist, totalscore, nfo);
//    }
//    store splice sites internally for online remapping
//    bl_storeSpliceEvent (space, seq, list, events, 0, 100, seqs, querydesc);

  } else { 
    
//    fprintf(stdout, "destructing list\n");
    bl_gmatchlistDestruct(space, list);
    list = bl_gmatchlistInit(space, maxedist);
    rmalign = 1;
  }

  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  
  for (i=0; i < r-l+1; i++) {
    FREEMEMORY(space, bd[i]);
  }
  FREEMEMORY(space, bd);

  for(i=0; i < cur->nooffragments; i++) {
    FREEMEMORY(space, lmv[i]);
    FREEMEMORY(space, lmr[i]);
    FREEMEMORY(space, lmc[i]);
    wrapMultiCharSeqAlignment(space, &a[i]);
  }
  
  FREEMEMORY(space, lmv);
  FREEMEMORY(space, lmr);
  FREEMEMORY(space, lmc);
  FREEMEMORY(space, aligns);
  FREEMEMORY(space, a);
  FREEMEMORY(space, b);
  FREEMEMORY(space, M);

  return list;
}


/*------------------------------ se_kdSplitRead ------------------------------
 *    
 * @brief find the splits of a chimeric reads from matchstem data
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdSplitRead(void *space, Suffixarray *arr, MultiCharSeq *seq, 
    char *querydesc, matchstem_t **stems, char **seqs, Uint len, 
    karlin_t *stats, spliceevents_t *events, segemehl_t *nfo) 
{
  int indel = -2;
  int scores[]={1, -2};
  Uint noofchains;
  gmatchlist_t* list;
  branchfragment_t* fragments;
  branchChain_t *chains;

  chains = branchChain(space, arr, stems, seqs, len, stats, 
      &noofchains, &fragments);

  list = se_kdAlignSplitChain (space, chains, noofchains,
      arr, seq, querydesc, stems, seqs, len, scores, indel, events, nfo);

  wrapChains(space, chains, noofchains);
  FREEMEMORY(space, fragments);
  FREEMEMORY(space, chains);

  return list;
}


/*--------------------------------- se_clip ----------------------------------
 *    
 * @brief clipping sequences
 * @author Steve Hoffmann 
 *   
 */
 
void
se_clip (void *space, fasta_t *reads, Uint elem, segemehl_t *nfo)
{

  if(nfo->hardclip3Prime || nfo->hardclip5Prime) {
    bl_fastaHardClip(space, reads, elem, nfo->hardclip5Prime, 
        nfo->hardclip3Prime);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateHardClip(space, reads, elem, nfo->hardclip5Prime, 
          nfo->hardclip3Prime);
    } 
  } 

  if(nfo->softclip3Prime || nfo->softclip5Prime) {
    bl_fastaSoftClip(space, reads, elem, 
        nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
        nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateSoftClip(space, reads, elem, 
          nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
          nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    }
  }

  return ;
}

/*----------------------------- se_kdGenomeMatch -----------------------------
 *    
 * @brief map reads to the genome
 * @author Steve Hoffmann 
 *   
 */

void
se_kdGenomeMatch(void *space, Suffixarray *s, fasta_t *reads, 
    segemehl_t *nfo) {

  unsigned char matchflag, matematchflag;
  matchstatus_t pairStatus = QUERY;
  char *seqs[2], *mateseqs[2];
  Uint k, i, u, *enctab, dim, wordno, len, matelen=0, jump, maxedist, 
       matemaxedist=0;
  karlin_t stats;
  bitvector *D, *Mv;
  Gmap map;
  gread_t read;
  matchstem_t *stems[2] = {NULL, NULL}, *matestems[2] = {NULL, NULL},
              *b0[2], *mateb0[2];
  gmatchlist_t *list=NULL, *matelist=NULL, *templist, 
               *bestpairlist=NULL, *slist=NULL, *slist2=NULL;
  spliceevents_t *events;
  bestseed_t best, bestmate;

  events = ALLOCMEMORY(space, NULL, spliceevents_t, 1);
  events->noofevents = 0;
  events->event = NULL;

  enctab = encodetab(nfo->seq->map, nfo->seq->mapsize);
  dim = reads->maxlen+1000;

  if(bl_fastaHasMate(reads)) {
    dim += nfo->maxinsertsize;
  }

  dim += 2*((reads->maxlen-ceil((nfo->accuracy*reads->maxlen)/100))+4);
  wordno = (reads->maxlen/BITVECTOR_WORDSIZE)+1; 

  D = ALLOCMEMORY(space, NULL, bitvector, 2*(dim+1));
  Mv = &D[dim+1];

  for(i=0; i <= dim; i++) {
    D[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mv[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  }  

  karlinunitcostpp(space, &stats.lambda, &stats.H, &stats.K);

  for (k=0; k < reads->noofseqs; k++) {
    pairStatus = QUERY;
    matchflag = 0;
    matematchflag = 0;
    best.mat = 0;
    bestmate.mat = 0;

    if(!nfo->mute) se_updateProgressBar(k, nfo);
    se_clip(space, reads, k, nfo);

    seqs[0] = bl_fastaGetSequence(reads, k);
    len = bl_fastaGetSequenceLength(reads, k);

#ifdef HASHING
    if (bl_fastaGetQuantity(reads, k) == 1){
      DBG("%u: %s\t%u\n",  k, bl_fastaGetSequence(reads, k), bl_fastaGetQuantity(reads, k));
    }

    //    fprintf(nfo->dev,"@%s\n%s\n+\n%s\n", 
    //        bl_fastaGetDescription(reads,k), seqs[0], bl_fastaGetQuality(reads,k));
    //    continue;
    continue;
#endif
    //    pthread_mutex_lock(nfo->mtx2);
    //    fprintf(nfo->dev, "%s\n", bl_fastaGetDescription(reads,k));
    //    fprintf(nfo->dev, "%s\n", bl_fastaGetMateDescription(reads,k));
    //    pthread_mutex_unlock(nfo->mtx2);


    if(len >= nfo->minsize) {  
      seqs[1] = charIUPACcomplement(space, seqs[0], len);

      /* convert for seed search */
      if (nfo->bisulfite){    
        seqs[0] = ALLOCMEMORY(space, NULL, char, len+1);
        memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
        bl_convertBisulfite(seqs, len, nfo->bisulfite, 1);
      }

      initGmap(&map, nfo->seq, 1);
      initRead(&read, k);

      if (nfo->jump == 0) {
        jump = floor(len/75) * 2;
        jump = (jump > 0) ? jump : 1;
      } else {
        jump = nfo->jump;
      }

      stems[0] = NULL; stems[1] = NULL;
      b0[0] = NULL; b0[1] = NULL;

      /* restrict search to one strand */
      for (u = 0; u < 2; u++){
        /* nfo->strand == 1 : search only on plus strand
         * => init stems[1] as empty
         * nfo->strand == 2 : search only on minus strand
         * => init stems[0] as empty
         * Note: empty but initialized stems are ignored
         * in function kdbest
         */
        if (nfo->strand == 2 - u){	 
          stems[u] = ALLOCMEMORY(space, NULL, matchstem_t, len);
          for (i = 0; i < len; i++){
            stems[u][i].branches = NULL;
            stems[u][i].noofbranches = 0;
          }
        }
      }

      /*
       * try to find full match, only possible
       * if there are no more than k_p
       * unmatchable characters are in read
       */
      if (nfo->bestonly && countNonMatchingChars(seqs[0], len) <= nfo->k_p){
        kdbest(space, s, seqs, len, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, stems, b0);
      }

      if (stems[0] == NULL){
        stems[0]=kdseeds(space, s, seqs[0], len, jump, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, b0[0]);
      }
      if (stems[1] == NULL){
        stems[1]=kdseeds(space, s, seqs[1], len, jump, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, b0[1]);
      }

      /* convert for alignment */
      if (nfo->bisulfite){
        FREEMEMORY(space, seqs[1]);
        memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
        seqs[1] = charIUPACcomplement(space, seqs[0], len);
        bl_convertBisulfite(seqs, len, nfo->bisulfite, 0); 
      }

      list = se_kdMatchStemAlign(space, s, nfo->seq, stems, seqs,
          len, &stats, nfo, enctab, D, &best);

      if (bl_fastaHasMate(reads)) {
        bestpairlist = NULL;

        mateseqs[0] = bl_fastaGetMate(reads, k);
        matelen = bl_fastaGetMateLength(reads, k); 

        mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen); 

        /* convert for direct mate alignment */
        if (nfo->bisulfite){
          mateseqs[0] = ALLOCMEMORY(space, NULL, char, matelen+1);
          memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
          bl_convertBisulfite(mateseqs, matelen, nfo->bisulfite, 0);
        }   

        if(se_kdMatchListhasMatches(list)) {
          matemaxedist = matelen-ceil((nfo->accuracy*matelen)/100);
          se_kdAlignMate(space, nfo->seq, mateseqs, matelen, 
              list, matemaxedist, enctab, D, nfo->maxinsertsize);
        }

        if (se_kdMatchListhasMatches(list) &&
            se_kdMatchListhasMates(list)) {
          /*pair is fully matched*/
          matematchflag = 1;
          pairStatus = PAIR;
        } else { 
          /*try to find mate first*/
          if (nfo->jump == 0) {
            jump = floor(matelen/75) * 2;
            jump = (jump > 0) ? jump : 1;
          } else {
            jump = nfo->jump;
          }

          /* convert for mate seed search */
          if (nfo->bisulfite){
            FREEMEMORY(space, mateseqs[1]);
            memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
            mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
            bl_convertBisulfite(mateseqs, matelen, nfo->bisulfite, 1);
          } 

          matestems[0] = NULL; matestems[1] = NULL;
          mateb0[0] = NULL; mateb0[1] = NULL;

          /* restrict search to one strand */
          for (u = 0; u < 2; u++){
            if (nfo->strand == 2 - u){
              matestems[u] = ALLOCMEMORY(space, NULL, matchstem_t, matelen);
              for (i = 0; i < matelen; i++){
                matestems[u][i].branches = NULL;
                matestems[u][i].noofbranches = 0;
              }
            }
          }

          /* 
           * try to find full match, only possible
           * if there are no more than k_p
           * unmatchable characters are in read
           */
          if (nfo->bestonly && countNonMatchingChars(mateseqs[0], matelen) <= nfo->k_p){
            kdbest(space, s, mateseqs, matelen, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, matestems, mateb0);
          }

          if (matestems[0] == NULL){
            matestems[0]=kdseeds(space, s, mateseqs[0], matelen, 
                jump, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, mateb0[0]);
          }
          if (matestems[1] == NULL){
            matestems[1]=kdseeds(space, s, mateseqs[1], matelen, 
                jump, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, mateb0[1]);
          }

          /* convert for mate alignment */
          if (nfo->bisulfite){
            FREEMEMORY(space, mateseqs[1]);
            memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
            mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
            bl_convertBisulfite(mateseqs, matelen, nfo->bisulfite, 0);
          }   

          matelist = se_kdMatchStemAlign(space, s, nfo->seq, matestems, 
              mateseqs, matelen, &stats, nfo, enctab, D, &bestmate);

          maxedist = len-ceil((nfo->accuracy*len)/100);

          se_kdAlignMate(space, nfo->seq, seqs, len, matelist, 
              maxedist, enctab, D, nfo->maxinsertsize);

          if (se_kdMatchListhasMatches(matelist) && 
              !se_kdMatchListhasMates(matelist) &&
              !se_kdMatchListhasMatches(list)) {
            /*query remains unmatched*/
            matematchflag = 1;
            pairStatus = MATE;
          }

          if(!se_kdMatchListhasMatches(matelist) &&
              se_kdMatchListhasMatches(list)) {
            /*mate remains unmatched*/
            pairStatus = QUERY;
          }

          if(se_kdMatchListhasMatches(list) &&
              se_kdMatchListhasMatches(matelist) && 
              !se_kdMatchListhasMates(matelist)) {
            /*pair not aligned properly but we have hits (long indel!)*/
            matematchflag = 1;
            bestpairlist = se_kdFindBestMatePair(space, list, matelist);
            pairStatus = PAIR_INS;
          } 

          if (se_kdMatchListhasMatches(matelist) && 
              se_kdMatchListhasMates(matelist)) {
            /*pair is fully matched in reverse order*/
            matematchflag = 1;
            templist = list;
            list = matelist;
            matelist = templist;
            pairStatus = PAIR_REV;
          }
        }
      }

      if (nfo->bestonly) {
        maxedist = list->minedist;
        if(matelist) matemaxedist = matelist->minedist;
      } else {
        maxedist = len-ceil((nfo->accuracy*len)/100);
        if(matelist) matemaxedist = matelen-ceil((nfo->accuracy*matelen)/100);
      }

      matchflag = 0;
      setReads(&map, &read, 1);

      /*report: single ends, fully matched pairs*/
      if(!bl_fastaHasMate(reads) || pairStatus == PAIR_REV || 
          pairStatus == PAIR) {
        if (list->n[0] || list->n[1]) matchflag = 1;
        se_setMatches(space, &read, list, maxedist, nfo);
        reportMatch(space, &map, reads, nfo, pairStatus, pairStatus == PAIR_REV);
        se_destructMatches(space, &read); 
      }

      /*report: spliced single ends */
      if(nfo->split && !bl_fastaHasMate(reads) && 
          !se_kdMatchListhasMatches(list)) {

        slist = se_kdSplitRead(space, s, nfo->seq, 
            bl_fastaGetDescription(reads, k), 
            stems, seqs, len, &stats, events, nfo);

        if (slist->n[0] || slist->n[1]) matchflag = 1;
        se_setMatches(space, &read, slist, maxedist, nfo);
        reportMatch(space, &map, reads, nfo, pairStatus, 0);
        se_destructMatches(space, &read); 
        bl_gmatchlistDestruct(space, slist);
      }

      /*report: bestpair from two separately calculated match lists*/
      if(pairStatus == PAIR_INS && bestpairlist) {
        se_setMatches(space, &read, bestpairlist, maxedist, nfo);
        matchflag = 1;
        matematchflag = 1;
        reportMatch(space, &map, reads, nfo, pairStatus, 0);
        se_destructMatches(space, &read); 
        bl_gmatchlistDestruct(space, bestpairlist);
      }

      /*report: spliced unmatched mate pairs*/
      if(bl_fastaHasMate(reads) && 
          (pairStatus == MATE || pairStatus == QUERY)) {

        if (nfo->split) {

          slist = NULL;
          slist2 = NULL;

          if(!se_kdMatchListhasMatches(list)) {
            slist = se_kdSplitRead(space, s,  nfo->seq,
                bl_fastaGetDescription(reads, k), 
                stems, seqs, len, &stats, events, nfo);

            if (slist->n[0] || slist->n[1]) matchflag = 1;
          } else {
            matchflag = 1;
          }

          if(!se_kdMatchListhasMatches(matelist)) {
            slist2 = se_kdSplitRead(space, s, nfo->seq,
                bl_fastaGetMateDescription(reads, k),
                matestems, mateseqs, matelen, &stats, events, nfo);

            if (slist2->n[0] || slist2->n[1]) matematchflag = 1;
          } else {
            matematchflag = 1;
          }
 

          if(slist && se_kdMatchListhasMatches(slist) 
              && (!slist2 || !se_kdMatchListhasMatches(slist2))) {

            /*spliced query full mate*/ 
            if (se_kdMatchListhasMatches(matelist)) {

              pairStatus = QUERY_SPL_FULL_MATE;

              se_setMatches(space, &read, matelist, matemaxedist, nfo);
              reportMatch(space, &map, reads, nfo, pairStatus, 1);
              se_destructMatches(space, &read); 
              /*spliced query no mate*/
            } else {
              pairStatus = QUERY_SPL_NO_MATE;
            }

            se_setMatches(space, &read, slist, maxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 0);        
            se_destructMatches(space, &read); 
          }

          if(slist2 && se_kdMatchListhasMatches(slist2) && 
              (!slist || !se_kdMatchListhasMatches(slist))) { 
            /*spliced mate full query*/ 
            if (se_kdMatchListhasMatches(list)) {

              pairStatus = MATE_SPL_FULL_QUERY;
              se_setMatches(space, &read, list, maxedist, nfo);
              reportMatch(space, &map, reads, nfo, pairStatus, 0);
              se_destructMatches(space, &read); 
              /*spliced query no mate*/
            } else {
              pairStatus = MATE_SPL_NO_QUERY;
            }

            se_setMatches(space, &read, slist2, matemaxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read); 
          }

          if(slist && se_kdMatchListhasMatches(slist) && slist2 && 
              se_kdMatchListhasMatches(slist2)) {
            /*both spliced*/
            pairStatus = PAIR_SPL;

            se_setMatches(space, &read, slist, maxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 0);
            se_destructMatches(space, &read); 

            se_setMatches(space, &read, slist2, maxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read); 
          }

          if(slist)  bl_gmatchlistDestruct(space, slist);
          if(slist2) bl_gmatchlistDestruct(space, slist2);

        } else {

          if(list && se_kdMatchListhasMatches(list)) {
            matchflag = 1;
            se_setMatches(space, &read, list, maxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 0);
            se_destructMatches(space, &read); 
          }

          if(matelist && se_kdMatchListhasMatches(matelist)) {
            matematchflag = 1;
            se_setMatches(space, &read, matelist, matemaxedist, nfo);
            reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read); 
          }
        }
      }

      bl_kdMatchstemDestruct(space, stems[0], len);
      bl_kdMatchstemDestruct(space, stems[1], len);
      if(matestems[0]) {
        bl_kdMatchstemDestruct(space, matestems[0], matelen);
        bl_kdMatchstemDestruct(space, matestems[1], matelen); 
        matestems[0] = NULL;
        matestems[1] = NULL;
      }

      bl_gmatchlistDestruct(space, list);
      if (nfo->bisulfite){
        FREEMEMORY(space, seqs[0]);
      }
      FREEMEMORY(space, seqs[1]);  

      if (bl_fastaHasMate(reads)) { 
        if (nfo->bisulfite){
          FREEMEMORY(space, mateseqs[0]);
        }
        FREEMEMORY(space, mateseqs[1]);
      }
      if (matelist) {
        bl_gmatchlistDestruct(space, matelist);
        matelist = NULL;
      }
    }
    
    bl_kdReportUnmatched(space, reads, k, matchflag, matematchflag, nfo);
  }

  wrapBitmatrix(space, D, 2*(dim+1));
  FREEMEMORY(space, D);
  FREEMEMORY(space, enctab);
  FREEMEMORY(space, events);
  return;
}


/*--------------------------- bl_kdReportUnmatched ---------------------------
 *    
 * @brief dump the unmatched sequences to a device
 * @author Steve Hoffmann 
 *   
 */

void
bl_kdReportUnmatched (void *space, fasta_t *reads, Uint k, 
    unsigned char matchflag, unsigned char matematchflag, segemehl_t *nfo)
{

  if(nfo->nomatchdev) { 

    if (!matchflag && !bl_fastaHasMate(reads)) { 
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);

      if (!bl_fastaHasQuality(reads)){
        fprintf(nfo->nomatchdev, ">%s\n%s\n", 
            bl_fastaGetDescription(reads, k), bl_fastaGetSequence(reads, k)); 
      } else {	   
        fprintf(nfo->nomatchdev, "@%s\n%s\n+%s\n%s\n",		
            bl_fastaGetDescription(reads, k), bl_fastaGetSequence(reads, k),
            bl_fastaGetDescription(reads, k), bl_fastaGetQuality(reads, k));
      }

      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }


    if ((!matchflag || !matematchflag) && bl_fastaHasMate(reads)) {
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);

      if(!matchflag) {
        if (!bl_fastaHasQuality(reads)){
          fprintf(nfo->nomatchdev, ">%s\n%s\n", 
              bl_fastaGetDescription(reads, k), bl_fastaGetSequence(reads, k)); 
        } else {	   
          fprintf(nfo->nomatchdev, "@%s\n%s\n+%s\n%s\n",		
              bl_fastaGetDescription(reads, k), bl_fastaGetSequence(reads, k),
              bl_fastaGetDescription(reads, k), bl_fastaGetQuality(reads, k));
        }
      }

      if(!matematchflag) {
        if (!bl_fastaHasQuality(reads)){
          fprintf(nfo->nomatchdev, ">%s\n%s\n", 
              bl_fastaGetMateDescription(reads, k), bl_fastaGetMate(reads, k)); 
        } else {	   
          fprintf(nfo->nomatchdev, "@%s\n%s\n+%s\n%s\n",		
              bl_fastaGetMateDescription(reads, k), bl_fastaGetMate(reads, k),
              bl_fastaGetMateDescription(reads, k), bl_fastaGetMateQuality(reads, k));
        }
      }

      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }

  }

  return ;
}

