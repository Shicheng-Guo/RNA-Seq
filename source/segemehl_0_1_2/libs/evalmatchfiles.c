
/*
 *  evalmatchfiles.c
 *  evalutation/statistics
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.10.2010 20:26:00 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "sort.h"
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "vtprogressbar.h"
#include "fileio.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "debug.h"
#include "evalmatchfiles.h"
#include "biofiles.h"
#include "splicesites.h"


/*-------------------------- bl_matchfileFrameStats --------------------------
 *    
 * @brief get descriptive statistics for a frame
 * @author Steve Hoffmann 
 *   
 */
 

matchfileFrameStats_t *
bl_matchfileFrameStats (void *space, matchfileFrame_t *frame) {
  Uint j, pos=0, ch, qu, rp, mc, rss, *c, *f, *n, 
       *D, D_ymax=0, D_ymax_1=0, D_xmax=0, noofstarts=0;
  matchfileFrameStats_t* stats;
  int **s;
  double *e, *r, *v, *y, *d, x=0;

  stats = ALLOCMEMORY(space, NULL, matchfileFrameStats_t, 1);
  v = ALLOCMEMORY(space, NULL, double, 256);
  d = ALLOCMEMORY(space, NULL, double, 256);
  e = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, double, 256);
  r = ALLOCMEMORY(space, NULL, double, 256);
  y = ALLOCMEMORY(space, NULL, double, 256);
  f = ALLOCMEMORY(space, NULL, Uint, 256);
  D = ALLOCMEMORY(space, NULL, Uint, MAX_STARTSITES);

  memset(v, 0, sizeof(double)*256);
  memset(c, 0, sizeof(double)*256);
  memset(r, 0, sizeof(double)*256);
  memset(y, 0, sizeof(double)*256);
  memset(d, 0, sizeof(double)*256);
  memset(f, 0, sizeof(Uint)  *256);
  memset(D, 0, sizeof(Uint)*MAX_STARTSITES);

  for(j=0; j < 256; j++) e[j] = log10(0);

  for(pos=0; pos < frame->width; pos++) {

    x += frame->cs[pos].len;
    s = ALLOCMEMORY(space, NULL, int*, 256);
    n = ALLOCMEMORY(space, NULL, Uint, 256);
    memset(s, 0, sizeof(int*)*256);
    memset(n, 0, sizeof(Uint)*256);
    f[(int)frame->cs[pos].ref]++;
    rss = frame->cs[pos].starts;
    noofstarts += rss;

    if(rss >= MAX_STARTSITES) {
      rss = MAX_STARTSITES-1;
      D_xmax = MAX_STARTSITES;
    } else {
      D_xmax = (rss > D_xmax) ? rss: D_xmax;
    }
    
    D[rss]++;
    if(D[rss] > D_ymax) {
      D_ymax = D[rss];
    }

    if(rss > 0 && D[rss] >D_ymax_1) {
      D_ymax_1 = D[rss];
    }

    for(j=0; j < frame->cs[pos].len; j++) {
      ch = frame->cs[pos].chars[j];
      qu = frame->cs[pos].quals[j];
      rp = frame->cs[pos].readpos[j];
      mc = frame->cs[pos].matchcnt[j];


      if(ch != frame->cs[pos].ref) {
        d[ch]++;
      }

      s[ch] = ALLOCMEMORY(space, s[ch], Uint, n[ch]+1);
      s[ch][n[ch]] = rp;
      e[ch] = log10add(e[ch],(qu/-10.));
      r[ch] += rp;
      y[ch] += mc;
      c[ch]++;
      n[ch]++;
    }

    for(j=0; j < 256; j++) {
      if(n[j]) v[j] += sqrt(var_int(s[j],n[j]));
      FREEMEMORY(space, s[j]);
    }
    
    FREEMEMORY(space, n);
    FREEMEMORY(space, s);
  }

  for(j=0; j < 256; j++) {
    if(c[j]) e[j] -= log10(c[j]);
    if(c[j]) r[j] /= c[j];
    if(c[j]) y[j] /= c[j];
    if(c[j]) v[j] /= f[j];
    if(c[j]) d[j] /= c[j];
  }


  stats->ntcnt = c;
  stats->mean_err = e;
  stats->mean_sde = v;
  stats->mean_pos = r;
  stats->mean_mul = y;
  stats->mean_dis = d;
  stats->char_frq = f;
  stats->mean_cov = x/frame->width;
  stats->dist_rss_xmax = D_xmax;
  stats->dist_rss_ymax = D_ymax;
  stats->dist_rss_ymax_1 = D_ymax_1;
  stats->dist_rss = D;
  stats->rss = noofstarts;

  return stats;
}



/*------------------------ bl_matchfileDestructFrame -------------------------
 *    
 * @brief wrap the frame
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructFrame(void *space, matchfileFrame_t *frame)
{
  
  bl_matchfileDestructCross(space, frame->cs, frame->width);
  FREEMEMORY(space, frame->cs);
  FREEMEMORY(space, frame);

  return ;
}

/*---------------------- bl_matchfileDestructFrameStats ----------------------
 *    
 * @brief destruct the frame statistics structure
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileDestructFrameStats(void *space, matchfileFrameStats_t *stats) {
  FREEMEMORY(space, stats->mean_err);
  FREEMEMORY(space, stats->mean_sde);
  FREEMEMORY(space, stats->mean_pos);
  FREEMEMORY(space, stats->mean_mul);
  FREEMEMORY(space, stats->mean_dis);
  FREEMEMORY(space, stats->dist_rss);
  FREEMEMORY(space, stats->char_frq);
  FREEMEMORY(space, stats->ntcnt);
}

/*---------------------- bl_matchfileGetRegionConsError ----------------------
 *    
 * @brief get the error e of a of a region of cross section based 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */

  double
bl_matchfileGetRegionConsError (matchfileFrame_t *frame, 
    Uint pos, Uint range)
{
  Uint i, j, l, r;
  double e=0, ef=0;

  assert(frame->width > range);

  if(range/2 >= pos) {
    l = 0;
    r = range;
  } else { 
    if(range/2 + pos >= frame->width-1) {
      l = frame->width-1-range;
      r = frame->width-1;
    } else { 
      l = pos - range/2;
      r = pos + range/2-1;
    }
  }


  for(i=l; i < r; i++) { 

    for(j=0, e=0; j < frame->cs[i].len; j++) {
      if(frame->cs[i].chars[j] != frame->cs[i].cons)
        e++;
    }

    if(frame->cs[i].len) { 
      e /= (double)frame->cs[i].len;
      ef += e;
    }
  }
  
  return ef/(double) range;
}

/*---------------------- bl_matchfileGetCrossConsError -----------------------
 *    
 * @brief get the error e of a cross section based on the consensus 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileGetCrossConsError (matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  double e=0;
 
  for(j=0; j < frame->cs[pos].len; j++) {
    if(frame->cs[pos].chars[j] != frame->cs[pos].cons)
      e++;
  }

  if(frame->cs[pos].len) { 
    e /= frame->cs[pos].len;
  }

  return e;
}

/*-------------------------- bl_matchfileGetNTError --------------------------
 *    
 * @brief return average quality score, ie error probability,
 * in a cross section for each nucleotide present in the cross section
 * @author Steve Hoffmann 
 *   
 */

double*
bl_matchfileGetNTError(matchfileFrame_t *frame, Uint pos) {
  Uint j;
  Uint *c;
  double *p;


  p = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < 256; j++) {
    p[j] = log10(0);
  }

  for(j=0; j < frame->cs[pos].len; j++) {
    p[(int)frame->cs[pos].chars[j]] = 
      log10add(p[(int)frame->cs[pos].chars[j]],
          (((double)frame->cs[pos].quals[j])/-10.)); 
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) p[j] -= log10(c[j]);
  }

  FREEMEMORY(space, c);
  return p;
}

/*------------------------- bl_matchfileGetNTReadPos -------------------------
 *    
 * @brief get for each NT average qual in a cross section
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTQual(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += frame->cs[pos].quals[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}


/*------------------------- bl_matchfileGetNTReadPos -------------------------
 *    
 * @brief get for each NT average read position in a cross section
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTReadPos(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += frame->cs[pos].readpos[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}


/*----------------------- bl_matchfileGetNTReadPosVar ------------------------
 *    
 * @brief get for each NT read position variance in a cross section pos
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTReadPosVar(matchfileCross_t *cs) {
  double *v, **r;
  int *c, j;

  r = ALLOCMEMORY(space, NULL, double*, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  v = ALLOCMEMORY(space, NULL, double, 256);
  memset(r, 0, sizeof(double*)*256);
  memset(c, 0, sizeof(Uint)*256);
  memset(v, 0, sizeof(double)*256);

  for(j=0; j < cs->len; j++) {

    r[(int)cs->chars[j]] =
      ALLOCMEMORY(space, 
          r[(int)cs->chars[j]], double, 
          c[(int)cs->chars[j]]+1);
 
    r[(int)cs->chars[j]][c[(int)cs->chars[j]]]
      = (double)cs->readpos[j] / (double)cs->readlen[j];  

    c[(int)cs->chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) { 
      v[j] = var(r[j],c[j]);
    }
    FREEMEMORY(space, r[j]);
  }

  FREEMEMORY(space, r);
  FREEMEMORY(space, c);
  return v;
}


/*----------------------- bl_matchfileGetNTRedundancy ------------------------
 *    
 * @brief get for each NT average redundancy of reads (multiple read hit count) 
 * in the cross section
 * @author Steve Hoffmann 
 *   
 */

double*
bl_matchfileGetNTRedundancy(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(double)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += (int)frame->cs[pos].matchcnt[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}




/*-------------------------- bl_matchfileGetNTEdist --------------------------
 *    
 * @brief for each cross section in the frame:
 * Get the average edist of reads in the cross section
 * @author Steve Hoffmann 
 *   
 */
 
double*
bl_matchfileGetNTEdist(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *e;

  e = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(e, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);


  for(j=0; j < frame->cs[pos].len; j++) {
    e[(int)frame->cs[pos].chars[j]] 
      += (int)frame->cs[pos].edist[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) e[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return e;
}


/*------------------------- bl_matchfileGetNTCounts --------------------------
 *    
 * @brief for each cross section in the frame: 
 * Get the count for all nucleotides
 * @author Steve Hoffmann 
 *   
 */

Uint*
bl_matchfileGetNTCounts(matchfileCross_t *cs) {
  Uint j, *cnt;

  cnt = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(cnt, 0, sizeof(Uint)*256);

  for(j=0; j < cs->len; j++) {
    cnt[(int)cs->chars[j]]++;
  }

  return cnt;
}


/*------------------------- bl_matchfileGetConsensus -------------------------
 *    
 * @brief calculate the consensus bases in a frame
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileGetConsensus(matchfileFrame_t *frame) {
  Uint i, j, *cnt, max;
  char gapalign = 0;

  for(i=0; i < frame->width; i++) {
    if (frame->ref) {
      frame->cs[i].ref = frame->ref[i];
    }
    frame->cs[i].cons = '^';
    cnt = bl_matchfileGetNTCounts(&frame->cs[i]);

    if(frame->cs[i].len){
      max = uarraymax(cnt, 255);
      frame->cs[i].cons = (char)max;
    }

    FREEMEMORY(space, cnt);

    gapalign = 0;
    if(frame->cs[i].noofdels > 1) { 
      for(j=0; j < frame->cs[i].noofdels; j++) {
        if(frame->cs[i].dels[j].len > 1) {
          gapalign = 1;
        }
      }
      if(gapalign) 
      bl_matchfileGapAlign(frame->cs[i].dels, frame->cs[i].noofdels);
    }
  }

  return;
}


/*-------------------------- bl_matchfileCrossPrint --------------------------
 *    
 * @brief print a cross section to stderr
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossPrint(void *space, matchfileFrame_t *frame) {
  Uint i, j, width, start;
  matchfileCross_t *cs;

  cs = frame->cs;
  start = frame->start;
  width = frame->width;
  for(i=0; i < width; i++) {
    fprintf(stderr, "%d: %d\t%d\t%d\t%s\t%s\t", i, cs[i].len, 
        cs[i].starts, cs[i].ends, cs[i].chars, 
        cs[i].quals);
    for(j=0; j < cs[i].len; j++) {
      fprintf(stderr, "%d,", cs[i].row[j]);
    }
    fprintf(stderr, "\t");

    for(j=0; j < cs[i].len; j++) {
      fprintf(stderr, "%d,", cs[i].edist[j]);
    }
    fprintf(stderr, "\n");
  }
}



/*-------------------------- bl_matchfileSampleCmp ---------------------------
 *    
 * @brief cmp chromosome positions sampled below 
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileSampleCmp (Uint elemA, Uint elemB, void *toSort, void *info)
{
  PairUint *samples;
  samples = (PairUint*) toSort;

  if(samples[elemA].a > samples[elemB].a) 
    return 1;
  if(samples[elemA].a < samples[elemB].a)
    return 2;

  if(samples[elemA].b > samples[elemB].b)
    return 1;
  if(samples[elemA].b < samples[elemB].b)
    return 2;
	
  return 0;
}


/*----------------------- bl_matchfileInitSampleStats ------------------------
 *    
 * @brief initialize sample stats
 * @author Steve Hoffmann 
 *   
 */
 
matchfileSampleStats_t*
bl_matchfileInitSampleStats (void *space, Uint maxsample, Uint maxcover, Uint mincover, 
    Uint areasize, double maxareae)
{

  matchfileSampleStats_t *stats;

  stats = ALLOCMEMORY(space, NULL, matchfileSampleStats_t, 1);
  stats->n = maxsample;
  stats->px = .0;
  stats->e = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->mincover = mincover;
  stats->maxcover = maxcover;
  stats->areasize = areasize;
  stats->maxareae = maxareae;
  stats->e_N = 0;
  stats->px = .0;
  stats->V_N = 0;
  stats->V_mu = .1;
  stats->V_sd = .1;
  stats->V_ll = 0;
  stats->Vx_N = 0;
  stats->Vx_mu = .1;
  stats->Vx_sd = .1;
  stats->Vx_ll = 0;
  stats->P = 0;
  stats->X = 0;
  stats->N = 0;
    
  stats->RR_N = 0;
  stats->RR = ALLOCMEMORY(space, NULL, Uint, 11);
  memset(stats->RR, 0, sizeof(Uint)*11);

  stats->MM_N = 0;
  stats->MM = ALLOCMEMORY(space, NULL, Uint, 51);
  memset(stats->MM, 0, sizeof(Uint)*51);
  
  stats->e_mu = ALLOCMEMORY(space, NULL, double, 2);
  stats->e_sd = ALLOCMEMORY(space, NULL, double, 2);

  stats->e_mu[0] = 0.1;
  stats->e_mu[1] = 0.6;
  stats->e_sd[0] = 0.1;
  stats->e_sd[1] = 0.6;
  stats->e_ll =.0;

  stats->S = ALLOCMEMORY(space, NULL, double, 6*6);
  stats->S_N = ALLOCMEMORY(space, NULL, Uint, 6);
  stats->Sx = ALLOCMEMORY(space, NULL, double, 6*6);
  stats->Sx_N = ALLOCMEMORY(space, NULL, double, 6);
  stats->R = ALLOCMEMORY(space, NULL, Uint, 100*255);
  stats->R_N = ALLOCMEMORY(space, NULL, Uint, 100*255);
  stats->V = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->Vx = ALLOCMEMORY(space, NULL, double, maxsample); 
  stats->RP = ALLOCMEMORY(space, NULL, Uint, 100);
  stats->RP_N = ALLOCMEMORY(space, NULL, Uint, 100);
  stats->RQ = ALLOCMEMORY(space, NULL, Uint, 255);
  stats->RQ_N = ALLOCMEMORY(space, NULL, Uint, 255);
  memset(stats->S, 0, sizeof(double)*(6*6));
  memset(stats->S_N, 0, sizeof(Uint)*6);
  memset(stats->Sx, 0, sizeof(double)*(6*6));
  memset(stats->Sx_N, 0, sizeof(Uint)*6);
  memset(stats->R, 0, sizeof(Uint)*(100*255));
  memset(stats->RP, 0, sizeof(Uint)*(100));
  memset(stats->RQ, 0, sizeof(Uint)*(255));
  memset(stats->R_N, 0, sizeof(Uint)*(100*255));
  memset(stats->RP_N, 0, sizeof(Uint)*(100));
  memset(stats->RQ_N, 0, sizeof(Uint)*(255));

  memset(stats->V, 0, sizeof(double)*maxsample);
  memset(stats->Vx, 0, sizeof(double)*maxsample);

  return stats;
}


/*--------------------- bl_matchfileDestructSampleStats ----------------------
 *    
 * @brief destruct sample stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructSampleStats (void *space, matchfileSampleStats_t *stats)
{
  
  FREEMEMORY(space, stats->e);
  FREEMEMORY(space, stats->e_mu);
  FREEMEMORY(space, stats->e_sd);
  FREEMEMORY(space, stats->S);
  FREEMEMORY(space, stats->S_N);
  FREEMEMORY(space, stats->Sx);
  FREEMEMORY(space, stats->Sx_N);
  FREEMEMORY(space, stats->R);
  FREEMEMORY(space, stats->R_N);
  FREEMEMORY(space, stats->V);
  FREEMEMORY(space, stats->Vx);
  FREEMEMORY(space, stats->RP);
  FREEMEMORY(space, stats->RQ);
  FREEMEMORY(space, stats->RP_N);
  FREEMEMORY(space, stats->RQ_N);
  FREEMEMORY(space, stats->RR);
  FREEMEMORY(space, stats->MM);

  return ;
}

/*----------------------- bl_matchfileGetConditionals ------------------------
 *    
 * @brief get conditionals
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileGetConditionals (void *space, matchfileFrame_t *frame,
    Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{
  Uint i, j, varcnt=0, *c;
  double *s, *ntreadposvar, avgvar =.0, er, e;
  char NT[] ={'A','C','G','T'};
  matchfileSampleStats_t *stats = (matchfileSampleStats_t *) nfo;
  int rpos;

  
  /*nucleotide noise matrix*/
  /*normalize readpos to 100*/
  for(j=0; j < frame->cs[pos].len; j++) { 

    rpos = trunc(((double)(((double)frame->cs[pos].readpos[j]*100.0)/((double)frame->cs[pos].readlen[j]))));
    
    MATRIX2D(stats->R_N, 255, 
        //(int)frame->cs[pos].readpos[j], 
        rpos,
        (int)frame->cs[pos].quals[j])++;
    
    stats->RP_N[rpos] += 1;
    
    //stats->RP[(int)frame->cs[pos].readpos[j]]++;
    stats->RQ_N[(int)frame->cs[pos].quals[j]]+=1;

    if(frame->cs[pos].chars[j] != frame->cs[pos].cons) { 
      MATRIX2D(stats->R, 255, 
          //(int)frame->cs[pos].readpos[j], 
          rpos,
          (int)frame->cs[pos].quals[j])++;
       
      stats->RP[rpos] += 1;

      //stats->RP[(int)frame->cs[pos].readpos[j]]++;
      stats->RQ[(int)frame->cs[pos].quals[j]]++;
    }  
  }

  if(frame->cs[pos].len < stats->mincover 
      || frame->cs[pos].len > stats->maxcover) {
    return;
  }
     
  er = bl_matchfileGetRegionConsError(frame, pos, stats->areasize);
  if(er > stats->maxareae) {
    return;
  }
  
  e = bl_matchfileGetCrossConsError(frame, pos);
 
  /*substitution matrix*/
  if(e-er > stats->px) {
    c = stats->Sx_N;
    s = stats->Sx;
    stats->X++;
  } else {
    c = stats->S_N;
    s = stats->S;
    stats->P++;
  }
  stats->N++;

  for(j=0; j < frame->cs[pos].len; j++) {
    MATRIX2D(s, 6, (int) ntcode[(int)frame->cs[pos].cons],  
        (int) ntcode[(int) frame->cs[pos].chars[j]])++;
    c[(int)ntcode[(int) frame->cs[pos].cons]] += 1;
    
    if(e-er < stats->px){
      
      if(frame->cs[pos].edist[j] <= 10)
        stats->RR[(int)frame->cs[pos].edist[j]]++;
      else
        stats->RR[10]++;
      stats->RR_N++;

      if(frame->cs[pos].matchcnt[j] <= 50)
        stats->MM[frame->cs[pos].matchcnt[j]]++;
      else
        stats->MM[50]++;
      stats->MM_N++;
    }
  }
 
  /*read pos variance density for variant and non-variant nts*/
  ntreadposvar = bl_matchfileGetNTReadPosVar(&frame->cs[pos]);

  for(i=0; i < 4; i++ ){
    if(NT[i] != frame->cs[pos].cons && ntreadposvar[(int)NT[i]] > 0) {
      varcnt++;
      avgvar += ntreadposvar[(int)NT[i]];
    }
  }

  if(varcnt) {
    avgvar /= varcnt;
  }

  /*this is one option to extend the model*/
  if(stats->Vx_N < stats->n) { 
    stats->Vx[stats->Vx_N++] = avgvar;
  }

  if(stats->V_N < stats->n) {
    stats->V[stats->V_N++] = ntreadposvar[(int)frame->cs[pos].cons];
  }


  FREEMEMORY(space, ntreadposvar);
  return ;
}

/*------------------------- bl_matchfileSampleStats --------------------------
 *    
 * @brief get sample statistics
 * @author Steve Hoffmann 
 *   
 */
  
void
bl_matchfileGetErrorDensity(void *space, matchfileFrame_t *frame, 
    Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{  
  double e, er;
  matchfileSampleStats_t *stats = (matchfileSampleStats_t*) nfo; 

  if (frame->cs[pos].len < stats->mincover || 
      frame->cs[pos].len > stats->maxcover) 
    return;
  
  er = bl_matchfileGetRegionConsError(frame, pos, stats->areasize);
  
  if (er > stats->maxareae) return;
  
  e = bl_matchfileGetCrossConsError(frame, pos);

  if(e > 0 && stats->e_N < stats->n) {    
    stats->e[stats->e_N++]=e-er;
  }

  return ;
}



/*--------------------------- bl_matchfileSmallMap ---------------------------
 *    
 * @brief get a small map to quickly find expressed/covered sites
 * @author Steve Hoffmann 
 *   
 */
 
char**
bl_matchfileSmallMap (void *space, matchfile_t* file, Uint **mapsize)
{

  FILE *fp = NULL;
  stringset_t *fields = NULL;
  char *buffer = NULL, ch, *curchrom = NULL, *filename, **map = NULL;
  Uint buffersize = 1024, len = 0, curstart = 0, 
       curend = 0, i, j, bin, lastbin=0, id, lastid=-1, *msz,
       noofseqs =0;

  matchfileindex_t *index;
  unsigned char header = 1;
  unsigned char gzip, fmt;
  struct gzidxfile *gzf = NULL;
  
  filename = file->filename;
  gzip = file->gzip;
  fmt = file->fmt;
  index = file->index;

  noofseqs = index->noofchroms;

  if (gzip) {
    //gzindex = bl_zranGetIndex(filename, &gzlen);
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, file->index->gzindex, 0);
  } else {
    fp = fopen(filename, "r");
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  map = ALLOCMEMORY(space, NULL, char*, noofseqs);
  memset(map, 0, sizeof(char*)*noofseqs);
  msz = ALLOCMEMORY(space, NULL, Uint, noofseqs);
  memset(msz, 0, sizeof(Uint)*noofseqs);

  while((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch == '\n' && len > 0) {

      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len] = '\0';
      header = (header) ? bl_matchfileIsHeader(buffer, len, fmt) : header;

      if (!header) { 
        fields = tokensToStringset(space, "\t", buffer, len);
        curstart = bl_matchfileGetStartPos(fields, fmt);
        curend = bl_matchfileGetEndPos(fields, fmt);
        curchrom = bl_matchfileGetChrom(fields, fmt);
        
        if(curend+1 > 0 && curchrom) {   

          for(id=0; id < index->noofchroms; id++) {
            if(strcmp(curchrom, index->chromnames[id]) == 0) {
              break;
            }

            for(j=0; j < strlen(index->chromnames[id]); j++) {
              if (isspace(index->chromnames[id][j])) break;
            }

            if(strlen(curchrom) <= j && 
                strncmp(curchrom, index->chromnames[id], j) == 0) {
              break;
            }    
          }

          if(id != lastid) {
            lastbin = -1;
            lastid = id;
          }

          if (id >= index->noofchroms) {
            DBGEXIT("reference sequence '%s' not found\n", curchrom);
          }
          
          for(i=curstart; i < curend; i++) {      
            bin = i/255;
            if (bin != lastbin) {
              if(bin >= msz[id]) {
                map[id] = ALLOCMEMORY(space, map[id], char, bin+1);
                memset(&map[id][msz[id]], 0, sizeof(char)*((bin+1)-msz[id]));
                msz[id] = bin;
              }
            }
            if(map[id][bin]+1 > 0) map[id][bin]++;
            lastbin = bin;
          }

        }

        destructStringset(space, fields);
      }

      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      len = 0;

    } else {
      if(ch != '\n') buffer[len++] = ch;
    }
  }

  if(gzip) {
    bl_destructgzidxfile(gzf);
    FREEMEMORY(space, gzf);
  }

  FREEMEMORY(space, buffer);
  fclose(fp);

  *mapsize = msz;
  return map;
}

/*------------------ bl_evalmatchfileSampleCrossSections -------------------
 *    
 * @brief sample and execute f on it
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_matchfileSampleCrossSections(void *space, 
    matchfile_t *file, fasta_t *set, Uint n, 
    void (*f)(void *, matchfileFrame_t*, Uint, 
      matchfileFrameStats_t *, void *), void *info)
{
  PairUint *samplepos;
  Uint i=0, r=0, j=0, k, *cumchrlen, 
       *order, prev=0, nchr, curchrom=0, curstart=0,
       *mapsize = NULL;
  matchfileFrame_t *frame = NULL;   
  matchfileFrameStats_t *stats = NULL;
  Uint maxcover = 20000;
  Uint setsize = 10000000;
  char **maps;


  //init random number generator
  srand((unsigned)(time(0))); 
  nchr = file->index->noofchroms;

  samplepos = ALLOCMEMORY(space, NULL, PairUint, n+1);
  memset(samplepos, 0, sizeof(PairUint)*n+1);
  cumchrlen = ALLOCMEMORY(space, NULL, Uint, nchr);
  memset(cumchrlen, 0, sizeof(Uint)*nchr);

  MSG("generating small map\n");
  //sum up the length of the references (i.e. chromosomes)
  maps = bl_matchfileSmallMap (space, file, &mapsize);
  MSG("map generated\n");

  cumchrlen[0] = file->index->matchend[0] - file->index->matchstart[0] + 1;
  for(i=1; i < nchr; i++) {
    cumchrlen[i] = (file->index->matchend[i] - 
        file->index->matchstart[i] + 1) + cumchrlen[i-1];
  }
  
  //randomize n positions across the genome and deterimine their 
  //chromosomes
  i = 0;
  j = 0;
 
  while(i < n && j < setsize) {
     k=0;
     
     samplepos[i].b = 
       (Uint)(((double)cumchrlen[nchr-1]*rand()/RAND_MAX+1))+1;
     
     while(samplepos[i].b > cumchrlen[k] && k+1 < nchr) k++;   
     samplepos[i].a = k;
     
 //    if(k > 0) {
 //      fprintf(stderr, "violation: i=%d, samplepos[i].b=%d, k=%d, nchr=%d, cumchrlen[k]=%d\n",
  //          i, samplepos[i].b, k, nchr, cumchrlen[k]);
  //     for(j=0; j < k; k++) {
  //       fprintf(stderr, "cumchrlen[%d]=%d\n", j, cumchrlen[k]);
  //     }

  //     exit(-1);
  //   }

     prev = (k == 0) ? 0 : cumchrlen[k-1];

     if(    maps[samplepos[i].a] 
         && mapsize[samplepos[i].a] > (samplepos[i].b - prev)/255 
         && maps[samplepos[i].a][(samplepos[i].b - prev)/255] > 10) {
       i++;
       r++;
        }

     j++;
  }

  NFO("\n selected %d positions for sampling %d %d\n", i, j, n);

 
  for(i=0; i < nchr; i++) {
    if(maps[i]) { 
      FREEMEMORY(space, maps[i]);
    }
  }

  FREEMEMORY(space, maps);
  FREEMEMORY(space, mapsize);
 
  if(j == setsize) {
    DBG("current sample size %d is below the minimum\n", r);
    FREEMEMORY(space, samplepos);
    FREEMEMORY(space, cumchrlen);
    return 0;
  }
  
  for(i=0; i < n; i++) {
    assert(samplepos[i].a < 1);
  }

  //sort
  order = quickSort(space, samplepos, n, bl_matchfileSampleCmp, NULL);   

  initProgressBarVT();
  for(i=0; i < n; i++) {
    assert(samplepos[order[i]].a < 1);
  }
  
  //evaluate
  //to increase speed a frame of size FRAMESIZE is loaded
  for(i=0; i < n; i++) {

  
     progressBarVT("positions sampled.", n, i, 25);
    //is position on a new chromosome or in a new frame?
    if(!frame || samplepos[order[i]].a != curchrom || 
        samplepos[order[i]].b-prev+1 >= frame->start+frame->width) {

      if(frame) { 
        bl_matchfileDestructFrame(space, frame);
        frame = NULL;
        //bl_matchfileDestructFrameStats(space, stats);
      }

      curchrom = samplepos[order[i]].a;
      curstart = samplepos[order[i]].b;
      prev = (samplepos[order[i]].a == 0) ? 0 : cumchrlen[samplepos[order[i]].a-1];

      frame = bl_matchfileGetFrame(space, file, 
          file->index->chromnames[samplepos[order[i]].a], 
          curstart-prev+1, 20000, set, maxcover, NULL); 
      bl_matchfileGetConsensus(frame);
     // stats = bl_matchfileFrameStats (space, frame);
    }
    
    f(space, frame, samplepos[order[i]].b-curstart, stats, info);   
  }

  NFO("\n %d positions sampled.\n", n);
    
  if(frame) { 
    bl_matchfileDestructFrame(space,frame);
    frame = NULL;
//    bl_matchfileDestructFrameStats(space, stats);
  }

  FREEMEMORY(space, order);
  FREEMEMORY(space, samplepos);
  FREEMEMORY(space, cumchrlen);
  return 1;
}

/*------------------------ bl_matchfileCrossStatsInit ------------------------
 *    
 * @brief initalize cross stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsInit (matchfileCrossStats_t *css, Uint len)
{

  css->len = len;
  css->var_s = ALLOCMEMORY(space, NULL, double, len);
  css->var_rt = ALLOCMEMORY(space, NULL, double, len);
  css->var_rq = ALLOCMEMORY(space, NULL, double, len);
  css->var_rr = ALLOCMEMORY(space, NULL, double, len);
  css->var_rv = ALLOCMEMORY(space, NULL, double, len);
  css->var_mm = ALLOCMEMORY(space, NULL, double, len);
  css->sub = ALLOCMEMORY(space, NULL, double, len);

  return ;
}

/*---------------------- bl_matchfileDestructCrossStats ----------------------
 *    
 * @brief destruct cross stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsDestruct (matchfileCrossStats_t *css)
{

  if (css->var_s) FREEMEMORY(space, css->var_s);
  if (css->var_rt) FREEMEMORY(space, css->var_rt);
  if (css->var_rq) FREEMEMORY(space, css->var_rq);
  if (css->var_rr) FREEMEMORY(space, css->var_rr);
  if (css->var_rv) FREEMEMORY(space, css->var_rv);
  if (css->var_mm) FREEMEMORY(space, css->var_mm);
  if (css->sub) FREEMEMORY(space, css->sub);

  return ;
}

/*------------------------- bl_matchfilePrintVariant -------------------------
 *    
 * @brief print variant
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsPrint (matchfileCrossStats_t *css, matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref)
{

  double *nv;
  char *ch, *rq;
  unsigned char *ed;
  uint32_t *rp, *mc, i, rpos, *rl;

 
  rp = cs->readpos;
  rq = cs->quals;
  ed = cs->edist;
  ch = cs->chars;
  nv = bl_matchfileGetNTReadPosVar(cs);
  mc = cs->matchcnt;
  rl = cs->readlen;

  fprintf(stderr, "----stats----\n");
  
  for(i=0; i < cs->len; i++) { 

    rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));
    
    fprintf(stderr, "P(%c -> %c) = %f\n", ref, ch[i], css->var_s[i]);
    fprintf(stderr, "1-RP(%d)=%f\n", (int)rpos, css->var_rt[i]);
    fprintf(stderr, "1-RQ(%d)=%f\n", (int)rq[i], css->var_rq[i]);
    fprintf(stderr, "(1-RQ + 1-RP)/2=%f\n", logadd(css->var_rt[i], css->var_rq[i]) - log(2));
    fprintf(stderr, "RR(%d)=%f\n", ed[i], css->var_rr[i]);
    fprintf(stderr, "RV(%f)=%f (mu: %f sd:%f)\n", nv[(int)ch[i]], css->var_rv[i],
        stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd);
    fprintf(stderr, "MM(%d)=%f (%d/%d )\n", mc[i], css->var_mm[i], 
        stats->MM[MIN(mc[i],10)], stats->MM_N);
    fprintf(stderr, "subtotal: %f\n", css->sub[i]);
  }

  return ;
}


/*------------------------- bl_matchfileTestVariant --------------------------
 *    
 * @brief test variant
 * @author Steve Hoffmann 
 *   
 */

double
bl_matchfileTestNonVariant ( matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref, 
    matchfileCrossStats_t *css)
{

  double p = .0, *nv, maxrp = -100, maxrq = -100, minrq = 100, minrp = 100;
  char *ch, *rq;
  unsigned char *ed;
  Uint len, i, rpos;
  uint32_t *rp, *mc, *rl;

  len = cs->len;
  ch = cs->chars;
  rp = cs->readpos;
  rq = cs->quals;
  mc = cs->matchcnt;
  ed = cs->edist;
  rl = cs->readlen;

  nv = bl_matchfileGetNTReadPosVar(cs);

  for(i=0; i < 100; i++) {
    if(stats->RP[i]) { 
      maxrp = (maxrp < (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : maxrp;
  
      minrp = (minrp > (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : minrp;
    }
  }

  for(i=0; i < 255; i++) {
    if(stats->RQ[i]) { 
      maxrq = (maxrq < (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0))) ? 
        (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0)) : maxrq;

      minrq = (minrq > (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0))) ? 
        (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0)) : minrq;
    }
  }

  bl_matchfileCrossStatsInit(css, len);
  
  for(i=0; i < len; i++) { 
    if(stats->S_N[(int)ntcode[(int)ref]] && ntcode[(int)ch[i]] < 5 ) {

      css->var_s[i] = MAX(MINSUBPROB, log(MATRIX2D(stats->S, 6, (int)ntcode[(int)ref], (int)ntcode[(int)ch[i]])+1) 
          - log(stats->S_N[(int)ntcode[(int)ref]]));

      rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));

      if((int)ch[i] == (int)ref) {

        css->var_rt[i] = log(minrp) - log(((double)stats->RP[rpos] + 1.0)
            /((double) stats->RP_N[rpos] + 1.0));

        css->var_rq[i] = log(minrq) - log(((double)stats->RQ[(int)rq[i]] + 1.0)
            /((double) stats->RQ_N[(int)rq[i]] + 1.0));

        css->var_rr[i] = log(((double)stats->RR[MIN(ed[i],10)] + 1.0)/
            ((double)stats->RR_N + 1.0));

        if(nv[(int)ch[i]] > EPSILON_NV) {
          css->var_rv[i] = MAX(-5.0, log(univarnormcdf(nv[(int)ch[i]], stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd)));
        } else {
          css->var_rv[i] = -5.0; 
        }
        
        css->var_mm[i] = .0;

        if( stats->MM_N > stats->MM[MIN(mc[i],10)])
        css->var_mm[i] = log(((double)stats->MM[MIN(mc[i],10)]+1.0)/
            ((double)stats->MM_N + 1.0));

      } else {

        css->var_rt[i] = log(((double)stats->RP[rpos] + 1.0)
            /((double) stats->RP_N[rpos] + 1.0)) - log(maxrp);

        css->var_rq[i] = log(((double)stats->RQ[(int)rq[i]] + 1.0)
            /((double)stats->RQ_N[(int)rq[i]] + 1.0)) - log(maxrq);

        css->var_rr[i] = log(1.0-((double)stats->RR[MIN(ed[i],10)] + 1.0)/
            ((double)stats->RR_N + 1.0));

        if(nv[(int)ch[i]] > EPSILON_NV) {
          css->var_rv[i] = MAX(-5.0, log(1.0-univarnormcdf(nv[(int)ch[i]], 
                stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd))); 
        } else {
          css->var_rv[i] = -0.001;
        }

        css->var_mm[i] = .0;

        if( stats->MM_N > stats->MM[MIN(mc[i],10)])
        css->var_mm[i] = MAX(log(0.001), 
            log(1.0-(((double)stats->MM[MIN(mc[i],10)]+1.0)/
                ((double)stats->MM_N + 1.0))));

      }

      p += css->var_s[i];
      p += logadd(css->var_rt[i], css->var_rq[i]) - log(2);
      p += css->var_rr[i];
      p += css->var_rv[i];
      p += css->var_mm[i];
      css->sub[i] = p;
    }   
  }


  FREEMEMORY(space, nv);
  return p;
}


/*------------------------ bl_matchfileTestNonVariant ------------------------
 *    
 * @brief test non variant
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileTestVariant (matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref, matchfileCrossStats_t *css)
{ 
  
  double px = .0, *nv, maxrp = -100, maxrq = -100, minrq = 100, minrp = 100;
  char *ch, *rq;
  unsigned char *ed;   
  Uint len, i, curerr, rpos;
  uint32_t *rp, *mc, *rl;
 
  len = cs->len;
  ch = cs->chars;
  rp = cs->readpos;
  rq = cs->quals;
  mc = cs->matchcnt;
  ed = cs->edist;
  nv = bl_matchfileGetNTReadPosVar(cs);
  rl = cs->readlen;
  
  for(i=0; i < 100; i++) {
    if(stats->RP[i]) { 
      maxrp = (maxrp < (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : maxrp;
  
      minrp = (minrp > (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        (((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : minrp;
    }
  }

  for(i=0; i < 255; i++) {
    if(stats->RQ[i]) { 
      maxrq = (maxrq < (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0))) ? 
        (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0)) : maxrq;

      minrq = (minrq > (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0))) ? 
        (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0)) : minrq;
    }
  }

  bl_matchfileCrossStatsInit(css, len);
  
  /*variant*/
  for(i=0; i < len; i++) { 
   
    if(stats->Sx_N[(int)ntcode[(int)ref]] && ntcode[(int)ch[i]] < 5) {

      css->var_s[i] = MAX(MINSUBPROB, (log(MATRIX2D(stats->Sx, 6, (int)ntcode[(int)ref], (int)ntcode[(int)ch[i]])+1) 
        - log(stats->Sx_N[(int)ntcode[(int)ref]])));

      rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));
            
      css->var_rt[i] = log(minrp) - log(((double)stats->RP[rpos]+1.0)/
          ((double)stats->RP_N[rpos]+1.0));

      css->var_rq[i] = log(minrq) - log(((double)stats->RQ[(int)rq[i]] + 1.0)/
            ((double)stats->RQ_N[(int)rq[i]] + 1.0));

      curerr = ((int)ch[i] != (int)ref && ed[i]);

      css->var_rr[i] = log(((double)stats->RR[MIN(ed[i]-curerr,10)]+1.0)/
          ((double)stats->RR_N + 1.0));
        
      css->var_mm[i] = .0;

      if( stats->MM_N > stats->MM[MIN(mc[i],10)])
      css->var_mm[i] = log(((double)stats->MM[MIN(mc[i],10)]+1.0)/
          ((double)stats->MM_N + 1.0));

      /*read variance*/ 
      if(nv[(int)ch[i]] > EPSILON_NV) {
        css->var_rv[i] = MAX(-5.0, log(univarnormcdf(nv[(int)ch[i]], 
              stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd))); 
      } else {
        css->var_rv[i] = -5.0;
      }
       
      px += css->var_s[i];     
      px += logadd(css->var_rt[i], css->var_rq[i]) - log(2);
      px += css->var_rr[i];
      px += css->var_mm[i];
      px += css->var_rv[i];
      css->sub[i] = px;  
    }   
  }
	
  FREEMEMORY(space, nv);
  return px;
}



/*----------------------------- bl_matchfileTest -----------------------------
 *    
 * @brief test
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileTest(void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, char ref, 
    matchfileSampleStats_t *stats, unsigned char show, void *nfo)
{
  matchfileCrossStats_t css; 
  double p_cons = .0, p_consx = .0, p_ref = .0, p_refx = .0, PX, P; 
  Uint *cnt;
  
  cs->s_cons = 1;
  cs->s_consx = 0;
  cs->s_ref = 1;
  cs->s_refx = 0;

  if(cs->len < stats->mincover 
      || cs->len > stats->maxcover) return 0;
    
  if(cs->len && !cs->cons){
    cnt = bl_matchfileGetNTCounts(cs);
    cs->cons = (char) uarraymax(cnt, 255);
    FREEMEMORY(space, cnt);
  }
    
  stats->MM[10] = (int)((double)stats->MM_N*0.0001)+1;
  stats->RR[10] = (int)((double)stats->RR_N*0.0001)+1;

  PX = log((double)stats->X/stats->N);
  P = log((double)stats->P/stats->N);

  p_cons  = bl_matchfileTestNonVariant (cs, stats, cs->cons, &css);
  if (show)  bl_matchfileCrossStatsPrint (&css, cs, stats, cs->cons);
  bl_matchfileCrossStatsDestruct (&css);

  p_consx = bl_matchfileTestVariant (cs, stats, cs->cons, &css);
  if (show)  bl_matchfileCrossStatsPrint (&css, cs, stats, cs->cons); 
  bl_matchfileCrossStatsDestruct (&css);
  
  p_ref  = bl_matchfileTestNonVariant (cs, stats, ref, &css);
  if (show) bl_matchfileCrossStatsPrint (&css, cs, stats, ref);
  bl_matchfileCrossStatsDestruct (&css);

  p_refx = bl_matchfileTestVariant (cs, stats, ref, &css);
  if (show)  bl_matchfileCrossStatsPrint (&css, cs, stats, ref);
  bl_matchfileCrossStatsDestruct (&css);
  
  cs->p_cons = P+p_cons;
  cs->s_cons =P+p_cons-logadd(PX+p_consx, P+p_cons);
  cs->p_consx = PX+p_consx;
  cs->s_consx = PX+p_consx-logadd(PX+p_consx, P+p_cons);

  cs->p_ref = P+p_ref;
  cs->s_ref =P+p_ref-logadd(PX+p_refx, P+p_ref);
  cs->p_refx = PX+p_refx;
  cs->s_refx = PX+p_refx-logadd(PX+p_refx, P+p_ref);


  return 0;

}


/*----------------------- bl_matchfilePrintTestResult ------------------------
 *    
 * @brief dump the test result
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileTestPrint (matchfileFrame_t *f, Uint p)
{

  char type;
  
  if(f->cs[p].s_consx > f->cs[p].s_cons && 
      f->cs[p].s_refx > f->cs[p].s_ref) { 
    type = 'B';
  } else if (f->cs[p].s_consx > f->cs[p].s_cons) { 
    type = 'C';
  } else{ 
    type = 'R';
  }

  printf("%s\t%d\t%c\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n", 
      f->chrname, f->start+p, type, 
      f->ref[p], f->cs[p].cons, 
      f->cs[p].chars, f->cs[p].s_cons, 
      f->cs[p].s_consx, f->cs[p].s_ref, 
      f->cs[p].s_refx);

  return ;
}

/*----------------------- bl_matchfileGroupTestsReset ------------------------
 *    
 * @brief reset group tests
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsReset (matchfileTestGroups_t *g)
{

  memset(g->chars, 0, sizeof(char*)*g->noofgroups);
  memset(g->len, 0, sizeof(Uint)*g->noofgroups);
  memset(g->s_cons, 0, sizeof(double)*g->noofgroups);
  memset(g->s_ref, 0, sizeof(double)*g->noofgroups);
  memset(g->s_consx, 0, sizeof(double)*g->noofgroups);
  memset(g->s_refx, 0, sizeof(double)*g->noofgroups);
  memset(g->type, 0, sizeof(char)*g->noofgroups);
	
  return ;
}

/*------------------------ bl_matchfileGroupTestsInit ------------------------
 *    
 * @brief init the group tests
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileTestGroupsInit (void *space, matchfileTestGroups_t *g, Uint n)
{
  g->noofgroups = n;
  g->s_cons = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->s_ref = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->s_consx = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->s_refx = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->len = ALLOCMEMORY(space, NULL, Uint, g->noofgroups);
  g->chars = ALLOCMEMORY(space, NULL, char*, g->noofgroups);
  g->type = ALLOCMEMORY(sapce, NULL, char, g->noofgroups);
  
  bl_matchfileTestGroupsReset(g);

  return ;
}


/*---------------------- bl_matchfileTestGroupsDestruct ----------------------
 *    
 * @brief destruct the groups
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsDestruct (void *space, matchfileTestGroups_t *g)
{
  Uint i;
  for(i=0; i < g->noofgroups; i++) { 
    FREEMEMORY(space, g->chars[i]);
  }
 
  FREEMEMORY(space, g->chars);
  FREEMEMORY(space, g->len);
  FREEMEMORY(space, g->type);
  FREEMEMORY(space, g->s_cons);
  FREEMEMORY(space, g->s_ref);
  FREEMEMORY(sapce, g->s_consx);
  FREEMEMORY(space, g->s_refx);


  return ;
}


/*------------------------ bl_matchfileGroupAddResult ------------------------
 *    
 * @brief add a test result to a group
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsAddResult (matchfileTestGroups_t *g, Uint no, 
    matchfileCross_t *cs)
{
  char *chars = cs->chars;
  Uint len = cs->len;
  double s_cons = cs->s_cons;
  double s_ref = cs->s_ref;
  double s_consx = cs->s_consx;
  double s_refx = cs->s_refx;

  g->s_cons[no] = logadd(g->s_cons[no], s_cons);
  g->s_ref[no] = logadd(g->s_ref[no], s_ref);
  g->s_consx[no] = logadd(g->s_consx[no], s_consx);
  g->s_refx[no] = logadd(g->s_refx[no], s_refx); 

  g->chars[no] = ALLOCMEMORY(space, g->chars[no], char, g->len[no]+len+1);
  memmove(&g->chars[no][g->len[no]], chars, len);
  g->len[no]+= len;
  g->chars[no][g->len[no]] = 0;

  return ;
}



/*----------------------- bl_matchfileTestGroupsPrint ------------------------
 *    
 * @brief print group tests
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsPrint (matchfileTestGroups_t *g, Uint no, 
    matchfileFrame_t *f, Uint p)
{
  char cons;
  Uint *cnt, i;

  cnt = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(cnt, 0, sizeof(Uint)*256);

  for(i=0; i < g->len[no]; i++) {
    cnt[(int) g->chars[no][i]]++;
  }

  cons = (char) uarraymax(cnt, 255);
  FREEMEMORY(space, cnt);

  if(g->s_consx[no] > g->s_cons[no] && 
      g->s_refx[no] > g->s_ref[no]) { 
    g->type[no] = 'B';
  } else if (g->s_consx[no] > g->s_cons[no]) { 
    g->type[no] = 'C';
  } else { 
    g->type[no] = 'R';
  }

  printf("%s\t%d\t%c\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n", 
      f->chrname, f->start+p, 
      g->type[no], f->ref[p], cons, g->chars[no], 
      g->s_cons[no], g->s_consx[no], 
      g->s_ref[no], g->s_refx[no]);
  
  return ;
}

/*---------------------- bl_matchfileEvalCrossSections -----------------------
 *    
 * @brief evaluate all cross sections
 * @author Steve Hoffmann 
 *   
 */


void
bl_matchfileEvalCrossSections (void *space,  matchfile_t **files, 
    int *groups, Uint nooffiles, fasta_t *set, 
    Uint (*f)(void *, Uint fidx, Uint cidx, Uint pos, matchfileCross_t*, char ref, 
    matchfileSampleStats_t *, unsigned char, void *), void *nfo)
{
  Uint i, j, k, n=0, nchr = 0, curchrom = 0, pos=0, maxgroupno=0;
  Uint ccnt;
  matchfileFrame_t **frames = NULL;
  matchfileTestGroups_t *g = NULL;
  //unsigned char exclusive = 0;
  Uint maxcover = 20000;

  for(j=0; j < nooffiles; j++) {
    nchr = MAX(nchr, files[j]->index->noofchroms);
  }

  if(groups) {    
      
    g = ALLOCMEMORY(space, NULL, matchfileTestGroups_t, 1);
    for(i=0; i < nooffiles; i++) {
      maxgroupno = (groups[i]  > maxgroupno) ? groups[i] : maxgroupno;
    }
    maxgroupno+=1;
  }

  frames = ALLOCMEMORY(space, NULL, matchfileFrame_t*, nooffiles);
  memset(frames, 0, sizeof(matchfileFrame_t*)*nooffiles);

  for(k=0; k < nchr; k++) {  

    for(j=0, n=0; j < nooffiles; j++) {  
      n = MAX(n, files[j]->index->matchend[k]);
     }

    //evaluate
    //to increase speed a frame of size FRAMESIZE is loaded

    for(j=0; j < nooffiles; j++) {

      if(files[j]->index->matchend[k]>0) { 
        for(i=1; i < n; i++) {

          if(groups) { 
            bl_matchfileTestGroupsInit(space, g, maxgroupno);
          }

          //is position on a new chromosome or in a new frame?
          if(!frames[j] || k != curchrom || 
              i >= frames[j]->start + frames[j]->width) {

            curchrom = k;
            if(frames[j]) { 
              bl_matchfileDestructFrame(space, frames[j]);
              frames[j] = NULL;
            }

            if(files[j]->index->matchend[k] < i ||
                files[j]->index->matchstart[k] > i) continue;

            frames[j] = bl_matchfileGetFrame(space, files[j], 
                files[j]->index->chromnames[k], i, 10000000, set, maxcover, NULL); 
          }

          pos = i - frames[j]->start;


          if(frames[j]->cs[pos].len)
          f(space, j, k, i, &frames[j]->cs[pos], frames[j]->ref[pos], 
              files[j]->index->stats, 0, nfo);   

          if(frames[j]->cs[pos].s_consx > frames[j]->cs[pos].s_cons || 
              frames[j]->cs[pos].s_refx > frames[j]->cs[pos].s_ref) {

            if(!groups) { 
              bl_matchfileTestPrint (frames[j], pos);
            } 

            ccnt++;
          }

          if(groups && frames[j]->cs[pos].len){ 
            bl_matchfileTestGroupsAddResult (g, groups[j], &frames[j]->cs[pos]);
          }
        }
        if(frames[j]) {
          bl_matchfileDestructFrame(space, frames[j]);
          frames[j] = NULL;
        }

        //    if(groups) { 

        /*iter the groups*/
        //      for(j=0; j < maxgroupno; j++) {
        //        if(g->s_consx[j] > g->s_cons[j] || 
        //            g->s_refx[j] > g->s_ref[j]) {
        //          exclusive = 1;

        //         for(u=0; u < maxgroupno; u++) {
        //           if(u != j && (g->s_consx[u] > g->s_cons[u] || 
        //                 g->s_refx[u] > g->s_ref[u])) {
        //             exclusive = 0;
        //           }
        //         }

        //         if(exclusive) { 
        //           bl_matchfileTestGroupsPrint (g, j, frames[j], pos);
        //         }
        //       }
        //     }

        //     bl_matchfileTestGroupsDestruct (space, g);
        //   }
      }
    }
  }

  for(j=0; j < nooffiles; j++){
    if(frames[j]) bl_matchfileDestructFrame(space,frames[j]);
  }

  if(groups) {
    FREEMEMORY(space, g);
  }

  FREEMEMORY(space, frames);
  return ;
}


/*-------------------------- bl_matchfileDumpStats ---------------------------
 *    
 * @brief dump the stats
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileDumpSampleStats (matchfileSampleStats_t *stats)
{

  Uint i;
  //Uint j;

  fprintf(stderr, "edensity:\n");

  fprintf(stderr, "mu={%f }, sd=%f, ll=%f\n", 
      stats->e_mu[0], stats->e_sd[0], stats->e_ll);

  /*
     fprintf(stderr, "noise matrix\n");

     for(i=30; i < 100; i++) {
     for(j=0; j < 100; j++) {
     if(MATRIX2D(stats->R_N, 255, j, i))
     fprintf(stderr, "%d %d %f\n", i, j, 
     (double)MATRIX2D(stats->R, 255, j, i)/MATRIX2D(stats->R_N, 255, j, i));
     }
     }
     */

  fprintf(stderr, "P=%d, X=%d, N=%d\n", stats->P, stats->X, stats->N);

  fprintf(stderr, "readpos\n");
  for(i=0; i < 100; i++) {
    if(stats->RP_N[i])
      fprintf(stderr, "%d %d %d %f\n", i, stats->RP[i], stats->RP_N[i], log((double)stats->RP[i]/(double)stats->RP_N[i]));
  }

  fprintf(stderr, "readqual\n");
  for(i=0; i < 255; i++) {
    if(stats->RQ_N[i])
      fprintf(stderr, "%d %d %d %f\n", i, stats->RQ[i], stats->RQ_N[i], (double)stats->RQ[i] / stats->RQ_N[i]);
  }

  fprintf(stderr, "readerror\n");
  for(i=0; i < 11; i++) {
    fprintf(stderr, "%d %d %d %f\n", i, stats->RR[i], stats->RR_N, (double)stats->RR[i] / stats->RR_N);
  }

  fprintf(stderr, "multiple matches\n");
  for(i=0; i < 50; i++) {
    fprintf(stderr, "%d %d %d %f\n", i, stats->MM[i], stats->MM_N, (double)stats->MM[i] / stats->MM_N);
  }

  fprintf(stderr, "readstartvar\n");
  for(i=0; i < stats->V_N; i++) {
    fprintf(stderr, "%f\n", stats->V[i]);
  }
 
  fprintf(stderr, "readstartvar gaussian model\n");

  fprintf(stderr, "mu=%f, sd=%f, ll=%f\n", 
      stats->V_mu, stats->V_sd, stats->V_ll);
  
  fprintf(stderr, "readstartvar X");

  for(i=0; i < stats->Vx_N; i++) {
    fprintf(stderr, "%f\n", stats->Vx[i]);
  }
 
  fprintf(stderr, "readstartvar X gaussian model\n");

  fprintf(stderr, "mu=%f, sd=%f, ll=%f\n", 
      stats->Vx_mu, stats->Vx_sd, stats->Vx_ll);

  
  return ;
}
