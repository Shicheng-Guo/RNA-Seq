
/*
 *  splicesites.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06/15/2011 10:56:53 PM CEST
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
#include "list.h"
#include "biofiles.h"
#include "splicesites.h"
#include "matepairs.h"

/*--------------------- bl_matchfileGetDistantSplitSites ---------------------
 *    
 * @brief get the acceptor or donor split site. with type 'N' the function
 *        will return both
 * @author Steve Hoffmann 
 *   
 */

distsplitsites_t*
bl_matchfileGetDistantSplitSites (void *space, matchfileCross_t *cs, Uint pos,
    Uint cidx, char type, Uint *noofsites, Uint *checkptr)
{

  Uint i, j, n, k= 0;
  distsplitsites_t *sites = NULL;
  n = 0;

  for(i=0; i < cs->noofsplits; i++) {
    if(cs->splits[i].edgetype == type || type == 'N') {
    
      /*site already seen*/
      for(j=0; j < k; j++) {
        if(sites[j].distpos == cs->splits[i].edge &&
           sites[j].distcidx == cs->splits[i].edgechridx) {
          break;
        }
      }

      /*new site*/
      if(j==k) { 
        sites = ALLOCMEMORY(space, sites, distsplitsites_t, k+1);
        sites[k].cs = cs;
        sites[k].noofsplits = 0;
        sites[k].distpos = cs->splits[i].edge;
        sites[k].distcidx = cs->splits[i].edgechridx;
        sites[k].acceptor = 0;
        sites[k].donor = 0; 
        sites[k].pos = pos;
        sites[k].cidx = cidx;
        sites[k].seen = 0;
        sites[k].transsplits = 0;
        k++;
      } 
      
      if(cs->splits[i].edgetype == 'A') { 
        sites[j].acceptor++;
      } else { 
        sites[j].donor++;
      }

      if(cs->splits[i].trans) {
        sites[j].transsplits++;
      }

      sites[j].noofsplits++;

    }
  }

  *checkptr = 0;
  *noofsites = k;
  return sites;
}





/*----------------------- bl_matchfileDistSplitSiteCmp -----------------------
 *    
 * @brief compare dist splitsites
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileDistSplitSiteCmp (Uint no, void *list, void *elem,  void *nfo)
{

  distsplitsites_t *a, *b;
  List *l;
  int i;

  l = (List*) list;
  i = bl_listGetCur(list, no);

  a = (distsplitsites_t*) l->data;
  b = (distsplitsites_t*) elem;

  if(a[i].distcidx < b->distcidx) return 2;
  if(a[i].distcidx > b->distcidx) return 1;
  if(a[i].distpos < b->distpos) return 2;
  if(a[i].distpos > b->distpos) return 1;

  return 0;
}



/*---------------------- bl_matchfileDestructSpliceSite ----------------------
 *    
 * @brief destructs a splice site
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructSpliceSite (void *space, splicesite_t *s)
{
  
  FREEMEMORY(space, s->splitsites);
  FREEMEMORY(space, s->noofsplits);
  FREEMEMORY(space, s->leftsiteidx);
  FREEMEMORY(space, s->leftedgeweight);
  FREEMEMORY(space, s->leftedgeacceptor);
  FREEMEMORY(space, s->leftedgedonor);
  FREEMEMORY(space, s->rightedgeacceptor);
  FREEMEMORY(space, s->rightedgedonor);
  FREEMEMORY(space, s->rightsiteidx);
  FREEMEMORY(space, s->rightedgeweight);
  FREEMEMORY(space, s->cs);

  return ;
}


/*---------------------- bl_matchfileDestructSpliceMap -----------------------
 *    
 * @brief destruct splice map 
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructSpliceMap (void *space, splicemap_t *sm)
{

  Uint i;

  for(i=0; i < sm->noofsplicesites; i++) {
    bl_matchfileDestructSpliceSite(space, &sm->map[i]);
  } 

  for(i=0; i < (sm->interval)*2+1; i++) {
    FREEMEMORY(space, sm->charhist[i]);
    FREEMEMORY(space, sm->charhistA[i]);
  }

  FREEMEMORY(space, sm->histogram);
  FREEMEMORY(space, sm->chrcnt);
  FREEMEMORY(space, sm->charhist);
  FREEMEMORY(space, sm->charhistA);
  FREEMEMORY(space, sm->map);

  return ;
}


/*---------------------- bl_matchfileCondenseSpliceMap -----------------------
 *    
 * @brief condense the split map to a splice map
 * @author Steve Hoffmann 
 *   
 */
 
splicemap_t*
bl_matchfileSpliceMap (void *space, splitmap_t *map, Uint interval, 
    Uint minsplitno)
{

  Uint i, j, h, sum=0, last=0, k=0;
  distsplitsites_t **distsites, *iter, cursplice;
  Uint *noofdistsites;
  splicesite_t *s = NULL;
  splicemap_t *sm;
  Uint totalsplits, left, cur;
  Uint atypes, dtypes;
  Uint pstrands, mstrands;
  

  Uint transsplits;
  matchfileCross_t **cs;
  Uint *splitsites;
  Uint *noofsplits;
  List distlist; 

  Uint idx, check;
  char *ref = NULL; 
   
  sm = ALLOCMEMORY(space, NULL, splicemap_t, 1);
  sm->histogram = ALLOCMEMORY(space, NULL, Uint, interval*2+1);
  sm->chrcnt = ALLOCMEMORY(space, NULL, Uint, interval*2+1);
  sm->charhist = ALLOCMEMORY(space, NULL, Uint*, interval*2+1);
  sm->charhistA = ALLOCMEMORY(space, NULL, Uint*, interval*2+1);

  memset(sm->histogram, 0, sizeof(Uint)*(interval*2+1));
  memset(sm->chrcnt, 0, sizeof(Uint)*(interval*2+1));


  for(i=0; i < interval*2+1; i++) {
    sm->charhist[i] = ALLOCMEMORY(space, NULL, Uint, 255);
    sm->charhistA[i] = ALLOCMEMORY(space, NULL, Uint, 255);
    memset(sm->charhist[i], 0, 255*sizeof(Uint));
    memset(sm->charhistA[i], 0, 255*sizeof(Uint));
  } 

  bl_listInit(&distlist, 1000, sizeof(distsplitsites_t));

  /*
   * iter all splits
   */

  for(i=0; i < map->noofsplits; i++) {
       
    /*
     *  enter splice site in case
     *      a.  first split seen just to set last
     *      b.  new chromosome
     *      c.  current split i \notin [last,last+interval] 
     */

    if (i==0 || map->cidx[i] != map->cidx[last] || 
        map->pos[i] > map->pos[last]+interval) {

      if(i > 0) { 

        /*
         *  iter all splits in interval [last,i[ 
         */

        splitsites = ALLOCMEMORY(space, NULL, Uint, i-last+1);
        noofsplits = ALLOCMEMORY(space, NULL, Uint, i-last+1);
        distsites = ALLOCMEMORY(space, NULL, distsplitsites_t*, i-last+1);
        noofdistsites = ALLOCMEMORY(space, NULL, Uint, i-last+1);

        cs = ALLOCMEMORY(space, NULL, matchfileCross_t**, i-last+1);
        mstrands = 0;
        pstrands = 0;
        atypes = 0;
        dtypes = 0;
        totalsplits = 0;
        transsplits = 0;

        for(j=last; j < i; j++) {
          cs[j-last] = &map->cs[j];
          splitsites[j-last] = map->pos[j];
          noofsplits[j-last] = map->cs[j].noofsplits;
          totalsplits += map->cs[j].noofsplits;

          /*
           *  store array of distant split sites in array of arrays
           *  and number in array
           */

          distsites[j-last] = bl_matchfileGetDistantSplitSites (space, &map->cs[j], 
              map->pos[j], map->cidx[j], 'N', &noofdistsites[j-last], &check);

          for(h=0; h < map->cs[j].noofsplits; h++) {

            if(map->cs[j].splits[h].strand == '-') {
              mstrands++; 
            } else {
              pstrands++; 
            }

            if(map->cs[j].splits[h].edgetype == 'D') {
              dtypes++;
            } else {
              atypes++;
            }

            if(map->cs[j].splits[h].trans) {
              transsplits++;
            }
          }
        } 
       
        /*
         *  proceed only if the splicesite qualifies
         */

        if(totalsplits >= minsplitno) {

          /*
           *  alloc new splice site
           */

          s = ALLOCMEMORY(space, s, splicesite_t, k+1);
          
          s[k].cidx = map->cidx[last];
          s[k].chromname = bl_fastaGetDescription(map->set, s[k].cidx);
          s[k].start = map->pos[last];
          s[k].end = map->pos[i-1];
          s[k].noofsplitsites = i-last+1;
          s[k].splitsites = splitsites;
          s[k].noofsplits = noofsplits;
          s[k].cs = cs;
          s[k].totalsplits = totalsplits;
          s[k].atypes = atypes;
          s[k].dtypes = dtypes;
          s[k].pstrands = pstrands;
          s[k].mstrands = mstrands;
          s[k].transsplits = transsplits;
          s[k].type = (s[k].atypes > s[k].dtypes) ? 'A' : 'D';
          s[k].strand = (s[k].mstrands > s[k].pstrands) ? '-' : '+';
          
          s[k].noofleftsites = 0;
          s[k].leftsiteidx = NULL;
          s[k].leftedgeweight = NULL;
          s[k].leftedgeacceptor = NULL;
          s[k].leftedgedonor = NULL;
          s[k].lefttranssplits = NULL;
          
          s[k].noofrightsites = 0;
          s[k].rightsiteidx = NULL;
          s[k].rightedgeweight = NULL;
          s[k].rightedgeacceptor = NULL;
          s[k].rightedgedonor = NULL;
          s[k].righttranssplits = NULL;
          s[k].matesupport = NULL;
          
          /*
           *   get the median split
           */
          
          for(sum=0, j=last; j < i; j++) {
            if(sum < s[k].totalsplits/2) {
              sum += s[k].noofsplits[j-last];
            }        
            if(sum >= s[k].totalsplits/2) break;
          }

          assert(j < i);
          s[k].median = map->pos[j];
 
          /*
           *   insert dist splits to sorted list
           *   only insert if distpos is still to come
           *   no order on cidx yet
           */

          for(j=last; j < i; j++) {
            for(h=0; h < noofdistsites[j-last]; h++) { 
              if(distsites[j-last][h].distcidx != map->cidx[last] ||
                  distsites[j-last][h].distpos >= map->pos[i]) { 
                
                distsites[j-last][h].splicesite = k;

                bl_listBinarySearchInsert(&distlist, &distsites[j-last][h], 
                  bl_matchfileDistSplitSiteCmp, NULL);
              }
            }
          }

          /*
           * find and link the donor split sites
           */

          cursplice.distpos = (s[k].median > (Uint)ceil((double)interval/2.0)) ?
            s[k].median - (Uint)ceil((double)interval/2.0) : 0;
          cursplice.distcidx = map->cidx[last];

          left = binarySearch_left(&distlist, distlist.numofelem, &cursplice, 
              bl_matchfileDistSplitSiteCmp, NULL);
          
          left = bl_listGetCur(&distlist, left);
          iter = (distsplitsites_t*) distlist.data;

          while(left != -1) {
    
            if(iter[left].distcidx == map->cidx[last] &&
               iter[left].distpos  <= cursplice.distpos+interval&&
               iter[left].distpos  >= cursplice.distpos) {
            

              for(h=0; h < s[k].noofleftsites; h++) {
                if(iter[left].splicesite == s[k].leftsiteidx[h])
                  break;
              }

              if(h == s[k].noofleftsites) {
                
                s[k].leftsiteidx = 
                  ALLOCMEMORY(space, s[k].leftsiteidx, Uint, 
                      s[k].noofleftsites+1);

                s[k].leftsiteidx[h] = iter[left].splicesite;
 
                s[k].leftedgeweight = 
                  ALLOCMEMORY(space, s[k].leftedgeweight, Uint, 
                      s[k].noofleftsites+1);
                
                s[k].leftedgeweight[h] = 0;

                s[k].leftedgeacceptor = 
                  ALLOCMEMORY(space, s[k].leftedgeacceptor, uint16_t, 
                      s[k].noofleftsites+1);

                s[k].leftedgeacceptor[h] = 0;
                 
                s[k].leftedgedonor = 
                  ALLOCMEMORY(space, s[k].leftedgedonor, uint16_t, 
                      s[k].noofleftsites+1);

                s[k].leftedgedonor[h] = 0;
               
                s[k].lefttranssplits =
                  ALLOCMEMORY(space, s[k].lefttranssplits, uint16_t, 
                      s[k].noofleftsites+1);

                s[k].lefttranssplits[h] = 0;

                s[k].matesupport = 
                  ALLOCMEMORY(space, s[k].matesupport, uint16_t,
                      s[k].noofleftsites+1);

                s[k].matesupport[h] = 0;

                s[k].noofleftsites++;
              }
              
              s[k].leftedgeweight[h] += iter[left].noofsplits;
              s[k].leftedgeacceptor[h] += iter[left].acceptor;
              s[k].leftedgedonor[h] += iter[left].donor; 
              s[k].lefttranssplits[h] += iter[left].transsplits;

              iter[left].seen = 1;
   
              cur = left;
              left = distlist.nodes[left].next;
              //elem = bl_listUnlink(&distlist, cur, NULL);           
              //FREEMEMORY(space, elem);

            } else {
              break;
            }        
          }

          /* now update the donor site */

          for(h=0; h < s[k].noofleftsites; h++) {
            idx = s[k].leftsiteidx[h];
             
            s[k].matesupport[h] = bl_matchfileSearchMateLink(space, 
                  &(map->matemap),  s[k].median, s[k].cidx, 
                  s[idx].median, s[idx].cidx); 

            s[idx].matesupport =
              ALLOCMEMORY(space, s[idx].matesupport, uint16_t,
                  s[idx].noofrightsites+1);

            s[idx].matesupport[s[idx].noofrightsites] = 
              s[k].matesupport[h];

            s[idx].rightsiteidx =
              ALLOCMEMORY(space, s[idx].rightsiteidx, Uint, 
                  s[idx].noofrightsites+1);

            s[idx].rightsiteidx[s[idx].noofrightsites] = k;

            s[idx].rightedgeweight =
              ALLOCMEMORY(space, s[idx].rightedgeweight, Uint, 
                  s[idx].noofrightsites+1);

            s[idx].rightedgeweight[s[idx].noofrightsites] 
              = s[k].leftedgeweight[h];
 
            s[idx].rightedgeacceptor = 
                  ALLOCMEMORY(space, s[idx].rightedgeacceptor, uint16_t, 
                      s[idx].noofrightsites+1);

            s[idx].rightedgeacceptor[s[idx].noofrightsites] 
              = s[k].leftedgedonor[h];
                 
            s[idx].rightedgedonor = 
                  ALLOCMEMORY(space, s[idx].rightedgedonor, uint16_t, 
                      s[idx].noofrightsites+1);

            s[idx].rightedgedonor[s[idx].noofrightsites]  
              = s[k].leftedgeacceptor[h];

            s[idx].righttranssplits = 
                  ALLOCMEMORY(space, s[idx].righttranssplits, uint16_t, 
                      s[idx].noofrightsites+1);
            
            s[idx].righttranssplits[s[idx].noofrightsites]  
              = s[k].lefttranssplits[h];

            s[idx].noofrightsites++;
          }

          idx = bl_fastxFindIDIdx(s[k].chromname, map->set);
          ref = bl_fastaGetSequence(map->set, idx);

          /*
           *  just histograms
           */

          if(s[k].type == 'D') {
            if(s[k].strand=='+') {
              sm->charhist[interval+1][(int)ref[s[k].median]]++;
            } else {
              sm->charhist[interval+1][(int)charComplementChar(ref[s[k].median])]++;
            }
            for(h=1; h < 11; h++) { 
              if(s[k].strand == '+') {
                sm->charhist[interval+1-h][(int)ref[s[k].median-h]]++;  
                sm->charhist[interval+1+h][(int)ref[s[k].median+h]]++;
              } else {
                sm->charhist[interval+1-h][(int)charComplementChar(ref[s[k].median+h])]++;
                sm->charhist[interval+1+h][(int)charComplementChar(ref[s[k].median-h])]++;
              }
            }
          }

 
          if(s[k].type == 'A') {
            if(s[k].strand=='+') {
              sm->charhistA[interval+1][(int)ref[s[k].median]]++;
            } else {
              sm->charhistA[interval+1][(int)charComplementChar(ref[s[k].median])]++;
            }
            for(h=1; h < 11; h++) { 
              if(s[k].strand == '+') {
                sm->charhistA[interval+1-h][(int)ref[s[k].median-h]]++;  
                sm->charhistA[interval+1+h][(int)ref[s[k].median+h]]++;
              } else {
                sm->charhistA[interval+1-h][(int)charComplementChar(ref[s[k].median+h])]++;
                sm->charhistA[interval+1+h][(int)charComplementChar(ref[s[k].median-h])]++;
              }
            }
          }


          sm->histogram[interval+1] += s[k].noofsplits[j-last];
          sum = s[k].noofsplits[j-last];
          for(h=last; h < i; h++) {
            if(h > j) {
              sum +=  s[k].noofsplits[h-last];
              sm->histogram[interval+1+(h-j)] += s[k].noofsplits[h-last];
//              chrcnt[interval+1+(h-j)]++;
//              charhist[interval+1+(h-j)][(int)ref[s[k].median]+(h-j)]++;
            }
            if(h < j) {

              sum +=  s[k].noofsplits[h-last];
              sm->histogram[interval+1-(j-h)] += s[k].noofsplits[h-last];
//              chrcnt[interval+1-(j-h)]++;
//              charhist[interval+1-(j-h)][(int)ref[s[k].median-(j-h)]]++;
            }
          }
                k++;
        } else {

          FREEMEMORY(space, splitsites);
          FREEMEMORY(space, noofsplits); 
          FREEMEMORY(space, cs);                
        }

        /*
         *  since distsites are stored in a list -> destruct array
         */
         
        for(j=last; j < i; j++) {
          if(distsites[j-last])
            FREEMEMORY(space, distsites[j-last]);  
        }
        FREEMEMORY(space, distsites);
        FREEMEMORY(space, noofdistsites);
      }

      last = i;
    }
  }


  cur = distlist.first;
  while(cur != -1) { 
    iter = (distsplitsites_t *)bl_listGetElem(&distlist, cur);
    cur = distlist.nodes[cur].next;
  }
  bl_listDestruct(&distlist, NULL);
  
  sm->map = s;
  sm->noofsplicesites = k;
  sm->interval = interval;
 
 
  return sm;
}




/*----------------------- bl_matchfileDestructSplitMap -----------------------
 *    
 * @brief destruct the split map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructSplitMap (void *space, splitmap_t *splitmap)
{

  if(splitmap->pos) FREEMEMORY(space, splitmap->pos);
  if(splitmap->cidx) FREEMEMORY(space, splitmap->cidx);
 
  if(splitmap->cs) { 
    bl_matchfileDestructCross(space, splitmap->cs, splitmap->noofsplits);
    FREEMEMORY(space, splitmap->cs);
  }

  splitmap->pos = NULL;
  splitmap->cidx = NULL;
  splitmap->cs = NULL;
  bl_matchfileDestructMateBinMap(space, &(splitmap->matemap));

  return ;

}


/*------------------------ bl_matchfileInitSplicemap -------------------------
 *    
 * @brief init the splicing map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileInitSplitMap (void *space, splitmap_t *splitmap, annotationtrack_t *bed, 
    matchfile_t **files, fasta_t *set)
{
  
  splitmap->bed = bed;
  splitmap->noofsplits = 0;
  splitmap->files = files;
  splitmap->pos = NULL;
  splitmap->cidx = NULL;
  splitmap->cs = NULL;
  splitmap->set = set;

  bl_matchfileInitMateBinMap (space, &(splitmap->matemap));

  return ;
}


/*------------------------ bl_matchfileAddSplitToMap -------------------------
 *    
 * @brief adding splice site to splice map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileAddSplitToMap (void *space, splitmap_t *map, matchfileCross_t *cs, 
    Uint cidx, Uint pos, char ref)
{
  Uint k;

  k = map->noofsplits;

  map->cs = ALLOCMEMORY(space, map->cs, matchfileCross_t, k+1);
  map->pos = ALLOCMEMORY(space, map->pos, Uint, k+1);
  map->cidx = ALLOCMEMORY(space, map->cidx, Uint, k+1);

  bl_matchfileCopyCross(space, &map->cs[k], cs);
  map->pos[k] = pos;
  map->cidx[k] = cidx;

  map->noofsplits++;
  
	
  return ;
}



/*---------------------------- bl_matchfileSplit ----------------------------
 *    
 * @brief check for splits
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileSplit (void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs,
    char ref, matchfileSampleStats_t *stats, unsigned char show, void *nfo)
{
  splitmap_t *map = (splitmap_t*) nfo; 
  char *chr;

  /*handle groups!*/

  if(map) {
    
    chr = bl_fastaGetDescription(map->set, cidx);
    bl_matchfileAddToMateBinMap (space, &(map->matemap), cs, cidx, pos, ref);

    if(cs->noofsplits > 0) { 
      bl_matchfileAddSplitToMap(space, map, cs, cidx, pos, ref);
    }
    return 1;
  } else {
    fprintf(stderr, "No map");
  }
  	
  return 0;
}

/*------------------------------- printsplits --------------------------------
 *    
 * @brief print splits
 * @author Steve Hoffmann 
 *   
 */
 
void
printsplits (void *space, char *chr, Uint pos, matchfileCross_t *cs, 
    splitmap_t *map)
{
  Uint i, len, j, k=0;
  char *chr2, *ptr, *res;
  matchfileSplit_t *split;
  Uint *positions;
  char **chroms;
  Uint *poscnt;

  if(cs->noofsplits > 0) {
    printf("%s\t%d\t%d\t", chr, pos, cs->noofsplits);

    chroms = ALLOCMEMORY(space, NULL, char*, cs->noofsplits);
    positions = ALLOCMEMORY(space, NULL, Uint, cs->noofsplits);
    poscnt = ALLOCMEMORY(space, NULL, Uint, cs->noofsplits);

    memset(chroms, 0, sizeof(char*)*cs->noofsplits);
    memset(poscnt, 0, sizeof(Uint)*cs->noofsplits);
    memset(positions, 0, sizeof(Uint)*cs->noofsplits);

    for(i=0; i < cs->noofsplits; i++) {
      split = &cs->splits[i];
     
      len = bl_fastaGetDescriptionLength(map->set, split->edgechridx);
      chr2 = ALLOCMEMORY(space, NULL, char, len+1);
      memmove(chr2, bl_fastaGetDescription(map->set, split->edgechridx), len);
      chr2[len] = 0;

      for(j=0; j < k; j++) {
        if(chroms[j] && !strcmp(chroms[j], chr2) && positions[j] == split->edge) {
          poscnt[j]++;
          break;
        }
      }
      
      if(j == k) {
        chroms[j] = chr2;
        positions[j] = split->edge;
        poscnt[j] = 1;
        k++;
      } else {
        FREEMEMORY(space, chr2); 
      } 
    }


    for(j=0; j < k; j++) { 
      res = chroms[j];
      ptr = strtok(chroms[j], " ");
      printf("%s:%d:%d\t", ptr, positions[j], poscnt[j]);
      FREEMEMORY(space, res);
    }
    
    printf("\n");
    FREEMEMORY(space, chroms);
    FREEMEMORY(space, poscnt);
    FREEMEMORY(space, positions);
  }
	
  return ;
}

/*------------------------------ printsplicebed ------------------------------
 *    
 * @brief print the splice bed file
 * @author Steve Hoffmann 
 *   
 */

void
printsplicebed(void *space, splicemap_t *sm, char* title, FILE *dev, FILE *transdev) {
  Uint i, j;
  char *type;
 
  fprintf(dev, "track name=splicesites_%s description=\"splice sites %s\" useScore=0\n", title, title);
  fprintf(transdev, "track name=transsplicesites_%s description=\"trans splice sites %s\" useScore=0\n", title, title);


  for(i=0; i < sm->noofsplicesites; i++) {
    for(j=0; j < sm->map[i].noofrightsites; j++) { 
      if(sm->map[sm->map[i].rightsiteidx[j]].cidx == sm->map[i].cidx &&  
          sm->map[sm->map[i].rightsiteidx[j]].median - sm->map[i].median < 200000) { 
        
        if(sm->map[i].righttranssplits[j]) {
         
          fprintf(transdev, "%s\t%d\t%d\t%s:%d:%d:%d\t%d\t%c\n", 
              sm->map[i].chromname, sm->map[i].median, 
              sm->map[sm->map[i].rightsiteidx[j]].median,
              "transsplits", 
              sm->map[i].rightedgeweight[j],
              sm->map[i].righttranssplits[j],
              sm->map[i].totalsplits, 
             sm->map[i].matesupport[j],
//             (int)(500+((float)sm->map[i].totalsplits/50.0)*500),
              sm->map[i].strand);


        } else { 

          fprintf(dev, "%s\t%d\t%d\t%s:%d:%d\t%d\t%c\n", 
              sm->map[i].chromname, sm->map[i].median, 
              sm->map[sm->map[i].rightsiteidx[j]].median,
              "splits", 
              sm->map[i].rightedgeweight[j],
              sm->map[i].totalsplits,
              sm->map[i].matesupport[j],
//              (int)(500+((float)sm->map[i].totalsplits/50.0)*500),
              sm->map[i].strand);
        }
        
      } else {

        if(sm->map[i].righttranssplits[j]) 
          type = "transdistsplit";
        else
          type = "distsplit";
        
        fprintf(dev, "%s\t%d\t%d\t%s:%s:%d:%d:%d:%d\t%d\t%c\n", 
            sm->map[i].chromname, sm->map[i].median,  
            sm->map[i].median, 
            type,             
            sm->map[sm->map[i].rightsiteidx[j]].chromname,
            sm->map[sm->map[i].rightsiteidx[j]].median,
            sm->map[i].rightedgeweight[j],  
            sm->map[i].righttranssplits[j],
            sm->map[i].totalsplits, 
            sm->map[i].matesupport[j],
//            (int)(500+((float)sm->map[i].totalsplits/50.0)*500),
            sm->map[i].strand);
      }
    }

    for(j=0; j < sm->map[i].noofleftsites; j++) { 
      if(sm->map[sm->map[i].leftsiteidx[j]].cidx != sm->map[i].cidx
          || sm->map[i].median - sm->map[sm->map[i].leftsiteidx[j]].median >= 200000) { 
        
        if(sm->map[i].lefttranssplits[j]) 
          type = "transdistsplit";
        else
          type = "distsplit";
 
        fprintf(dev, "%s\t%d\t%d\t%s:%s:%d:%d:%d:%d\t%d\t%c\n", 
            sm->map[i].chromname, sm->map[i].median, sm->map[i].median,
            type, 
            sm->map[sm->map[i].leftsiteidx[j]].chromname,
            sm->map[sm->map[i].leftsiteidx[j]].median, 
            sm->map[i].leftedgeweight[j],
            sm->map[i].lefttranssplits[j],
            sm->map[i].totalsplits, 
            sm->map[i].matesupport[j],
 //           (int)(500+((float)sm->map[i].totalsplits/50.0)*500), 
            sm->map[i].strand);
      }
    }

  }
}


/*----------------------- bl_matchfileSpliceAnnotation -----------------------
 *    
 * @brief intersect splicesites with exon annotation
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileSpliceAnnotation (void *space, splicemap_t* sm, annotationtrack_t *track)
{

  Uint i, j, left, len;
  int ldist, rdist;
  annotationitem_t cur;
  char *attr;

  for(i=0; i < sm->noofsplicesites; i++) {

    cur.chromname = sm->map[i].chromname;
    cur.strand = sm->map[i].strand;
    cur.start = (sm->map[i].median < sm->interval) ? 0 : 
      sm->map[i].median - sm->interval;
    cur.end = sm->map[i].median + sm->interval;

    left = binarySearch_left(track, track->noofitems, &cur, 
        bl_annotationitem_cmp_track, NULL);
    /*
       fprintf(stdout, "searching:%s,%d(%d, %d[i:%d]) -> leftmost %s,%d,%d\n",
       cur.chromname, cur.start, sm->map[i].median, sm->map[i].totalsplits, i,
       track->items[left].chromname,
       track->items[left].start,
       track->items[left].end);
       */
    while(track->items[left].start <= cur.end) {

      if(sm->map[i].noofleftsites || sm->map[i].noofrightsites) { 

        ldist = (int)sm->map[i].median - track->items[left].start;
        rdist = (int)sm->map[i].median - track->items[left].end;


        if((abs(ldist) <= sm->interval || abs(rdist) <= sm->interval) &&
            (sm->map[i].strand == track->items[left].strand)) {


          len = snprintf(NULL, 0, "%s \"%c\"", "alignedsplicesite_type", 
              sm->map[i].type);
          attr = ALLOCMEMORY(space, NULL, char, len+1);
          snprintf(attr, len+1, "%s \"%c\"", "alignedsplicesite_type", 
              sm->map[i].type);
          bl_GFFAddAttribute(space, &track->items[left], attr, len);
          FREEMEMORY(space, attr);


          len = snprintf(NULL, 0, "%s %d", "alignedsplicesite_pos", sm->map[i].median);
          attr = ALLOCMEMORY(space, NULL, char, len+1); 
          snprintf(attr, len+1, "%s %d", "alignedsplicesite_pos", sm->map[i].median);
          bl_GFFAddAttribute(space, &track->items[left], attr, len);
          FREEMEMORY(space, attr);

          len = snprintf(NULL, 0, "%s %d", "alignedsplicesite_dist", 
              (abs(ldist) < abs(rdist)) ? ldist : rdist);
          attr = ALLOCMEMORY(space, NULL, char, len+1);
          snprintf(attr, len+1, "%s %d", "alignedsplicesite_dist", 
              (abs(ldist) < abs(rdist)) ? ldist : rdist);
          bl_GFFAddAttribute(space, &track->items[left], attr, len);
          FREEMEMORY(space, attr);


          for(j=0; j < sm->map[i].noofrightsites; j++) {
            len = snprintf(NULL, 0, "%s \"%s\" %d %d", "distsplicesite_pos", 
                sm->map[sm->map[i].rightsiteidx[j]].chromname,
                sm->map[sm->map[i].rightsiteidx[j]].median,
                sm->map[i].rightedgeweight[j]);
            attr = ALLOCMEMORY(space, NULL, char, len+1);

            snprintf(attr, len+1, "%s \"%s\" %d %d", "distsplicesite_pos", 
                sm->map[sm->map[i].rightsiteidx[j]].chromname,
                sm->map[sm->map[i].rightsiteidx[j]].median,
                sm->map[i].rightedgeweight[j]);

            bl_GFFAddAttribute(space, &track->items[left], attr, len);
            FREEMEMORY(space, attr);
          }

          for(j=0; j < sm->map[i].noofleftsites; j++) {
            len = snprintf(NULL, 0, "%s \"%s\" %d %d", "distsplicesite_pos", 
                sm->map[sm->map[i].leftsiteidx[j]].chromname,
                sm->map[sm->map[i].leftsiteidx[j]].median,
                sm->map[i].leftedgeweight[j]);
            attr = ALLOCMEMORY(space, NULL, char, len+1);

            snprintf(attr, len+1, "%s \"%s\" %d %d", "distsplicesite_pos", 
                sm->map[sm->map[i].leftsiteidx[j]].chromname,
                sm->map[sm->map[i].leftsiteidx[j]].median,
                sm->map[i].leftedgeweight[j]);

            bl_GFFAddAttribute(space, &track->items[left], attr, len);
            FREEMEMORY(space, attr);
          }

        } 
      }
      left++;
    }
  }


  return ;
}

/*------------------------------- printsplice --------------------------------
 *    
 * @brief print splice
 * @author Steve Hoffmann 
 *   
 */
 
void
printsplice (void *space, splicemap_t *sm, FILE *dev)
{
   
 Uint i, j, noofdistsites, check;
 distsplitsites_t *distsites;

 for(i=0; i < sm->noofsplicesites; i++) {
    fprintf(dev, "%d\t%d\t%d\t%d\t%d\t", 
        sm->map[i].cidx, sm->map[i].median, 
        sm->map[i].totalsplits, 
        sm->map[i].start, sm->map[i].end
        );

    for(j=0; j < sm->map[i].noofleftsites; j++) {
      fprintf(dev, "left:%d,%d,%d,(acceptor:%d, donor:%d)\t",
          sm->map[sm->map[i].leftsiteidx[j]].cidx,
          sm->map[sm->map[i].leftsiteidx[j]].median,
          sm->map[i].leftedgeweight[j],
          sm->map[i].leftedgeacceptor[j],
          sm->map[i].leftedgedonor[j]);
    }

    for(j=0; j < sm->map[i].noofrightsites; j++) {
      fprintf(dev, "right:%d,%d,%d,(acceptor:%d, donor:%d)\t",
          sm->map[sm->map[i].rightsiteidx[j]].cidx,
          sm->map[sm->map[i].rightsiteidx[j]].median,
          sm->map[i].rightedgeweight[j],
          sm->map[i].rightedgeacceptor[j],
          sm->map[i].rightedgedonor[j]);
    }

    if(sm->map[i].noofrightsites == 0  && sm->map[i].noofleftsites == 0) {
       distsites = bl_matchfileGetDistantSplitSites (space, *sm->map[i].cs, 0, 0, 
        'N', &noofdistsites, &check);

      for(j=0; j < noofdistsites; j++) {
       fprintf(dev, "unconfirmed:%d,%d,%d,(acceptor:%d, donor:%d)\t",  
           distsites[j].distcidx, 
           distsites[j].distpos, 
           distsites[j].noofsplits,
           distsites[j].acceptor, 
           distsites[j].donor);
      }
    }

    fprintf(dev, "\n");
  }
 
  fprintf(dev, "Precision Histogram\n");

  for(i=0; i < sm->interval*2+1; i++) {
    fprintf(dev, "%u\t%u\n", i, sm->histogram[i]);
  }
  fprintf(dev, "Char Histogram Donor\n");

  for(i=0; i < sm->interval*2+1; i++) {
    fprintf(dev, "%d\t%d\t%d\t%d\t%d\t%d\n", i, sm->chrcnt[i], 
        sm->charhist[i]['A'],
        sm->charhist[i]['C'],
        sm->charhist[i]['G'],
        sm->charhist[i]['T']);
    FREEMEMORY(space, sm->charhist[i]);
  }

  fprintf(dev, "Char Histogram Acceptor\n");
  for(i=0; i < sm->interval*2+1; i++) {
    fprintf(dev, "%d\t%d\t%d\t%d\t%d\t%d\n", i, sm->chrcnt[i], 
        sm->charhistA[i]['A'],
        sm->charhistA[i]['C'],
        sm->charhistA[i]['G'],
        sm->charhistA[i]['T']);
    FREEMEMORY(space, sm->charhistA[i]);
  }

  return ;
}



