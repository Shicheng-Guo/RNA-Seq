
/*
 *  genfasta.c
 *  generate random fasta files
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07.01.2010 19:18:07 CET
 *  
 */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "biofiles.h"
#include "charsequence.h"
#include "randseqs.h"
#include "manopt.h"
#include "mathematics.h"

unsigned char mute=0;

int
main(int argc, char **argv) {

  manopt_optionset optset;
  manopt_arg *unflagged;
  manopt_intconstraint biseqconstraint;
  manopt_dblconstraint rateconstraint;
  annotationtrack_t* mytrack, *transtrack;
  fasta_t *fasta = NULL;
  Uint *space = NULL;
  Uint reflen = 1000000;
  Uint minreadlen = 75;
  Uint maxreadlen = 100;
  Uint minqual = 33;
  Uint maxqual = 73;
  Uint mindist = 150;
  Uint maxdist = 200;
  Uint polyAlen = 0;
  Uint alphabetlen = 4;
  Uint fiveprimelen = 0;
  Uint threeprimelen = 0;
  Uint n = 100;
  Uint i;
  geneset_t *geneset = NULL, *transset = NULL;
  double acc = 0.95;
  double Pmis = 0.8;//CHANGED from 0.6!!!
  double Pins = 0.1;//CHANGED from 0.2!!!
  double Pdel = 0.1;//CHANGED from 0.2!!!
  Uint biseq = 0;
  double rate = 0.5;
  char *alphabet = "ACTG";
  char *sequence = NULL;
  char *revcomp = NULL;
  char *subjectfilename = NULL;
  char *readsfilename = NULL;
  char *matefilename = NULL;
  char *fiveprime = NULL;
  char *threeprime = NULL;
  unsigned char fastq=0, split=0, splice=0;
  FILE *subjectdev = NULL;
  FILE *readsdev = stdout;
  FILE *matedev = NULL;
  FILE *dev = stdout;
  Uint maxchildren = 10;
  Uint mincisdist = 100;
  Uint maxcisdist = 1000;
  double Pcis = 0.9;
  double Pstrandswitch = 0.5; 
  char *dbfilename = NULL;
  char *splicefilename = NULL;
  char spliceedges =0;
  Uint trans = 0;
  char simulate = 0;
  char isoforms = 0;
  char distsplice = 0;
  char transtype = 'S'; 
  Uint cov = 10;

  rateconstraint.min = 0;
  rateconstraint.max = 1;
  biseqconstraint.min = 0;
  biseqconstraint.max = 2;

  srand((unsigned int)time(NULL));

  manopt_initoptionset(&optset, argv[0], NULL, 
      "  Generate random fasta sequences and simulated reads\n",
      "  GENFASTA is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      "0.1" ,
      "  Please report bugs to steve@bioinf.uni-leipzig.de");

  manopt(&optset, REQSTRINGOPT, 0, 'f', "readfile", 
      "path/filename to write the output to", "<file>", NULL, &readsfilename);
  manopt(&optset, REQSTRINGOPT, 0, 's', "subjectfile", 
      "path/filename to write subject sequence to", "<file>", NULL, &subjectfilename);
  manopt(&optset, FLAG, 0, 't', "split", 
      "split/spliced reads", NULL, NULL, &split);
  manopt(&optset, FLAG, 0, 'S', "splice", "generate a set of splice sites",
      NULL, NULL, &splice);
  manopt(&optset, REQSTRINGOPT, 0, 'p', "pairs", 
      "path/filename to write mate pair sequences to", NULL, NULL, &matefilename);
  manopt(&optset, REQUINTOPT, 0, 'l', "minreadlen", 
      "minimum size of queries", "<n>", NULL, &minreadlen);
  manopt(&optset, REQUINTOPT, 0, 'm', "maxreadlen", 
      "maximum size of queries", "<n>", NULL, &maxreadlen);
  manopt(&optset, FLAG, 0, 'q', "fastq", 
      "generate fastq reads", NULL, NULL, &fastq);
  manopt(&optset, REQSTRINGOPT, 0, 'a', "alphabet", 
      "alphabet for the fasta sequences", "<string>", NULL, &alphabet);
  manopt(&optset, REQINTOPT, 0, 'n', "readnumber", 
      "number of reads to be generated", "<string>", NULL, &n);
  manopt(&optset, REQINTOPT, 0, 'r', "reflen", 
      "length of reference sequence", "<string>", NULL, &reflen);
  manopt(&optset, DBLOPT, 0, 'A', "accuracy", 
      "accuracy of reads (1/errorrate)","<float>", NULL, &acc);
  manopt(&optset, DBLOPT, 0, 'M', "mismatches",
      "probablity of a mismatch for an erroneous site", "<float>", NULL, &Pmis);
  manopt(&optset, DBLOPT, 0, 'I', "insertions",
      "probablity of insertion for an erroneous site that is not a mismatch", "<float>", NULL, &Pins);
  
  manopt(&optset, REQUINTOPT, 0, 'B', "biseq",
      "bisulfite sequencing protocol (0 = no bisulfite, 1 = Lister et al., 2 = Cokus et al.)",
	 "<n>", &biseqconstraint, &biseq);
  manopt(&optset, DBLOPT, 0, 'R', "biseqrate", 
      "bisuflite conversion rate in reads","<float>", &rateconstraint, &rate);
  manopt(&optset, REQSTRINGOPT, 0, '5', "5prime",
      "add 5' adapter", "<string>", NULL, &fiveprime);
  manopt(&optset, REQSTRINGOPT, 0, '3', "3prime",
      "add 3' adapter", "<string>", NULL, &threeprime);
  manopt(&optset, REQINTOPT, 0, 'T', "polyA",
      "attach polyA tail", "<n>", NULL, &polyAlen);
  
  manopt(&optset, REQSTRINGOPT, 0, 'd', "database", 
      "path/filename of database sequences", "<file>", NULL, &dbfilename);
  manopt(&optset, REQSTRINGOPT, 0, 'b', "splicesites", 
      "path/filename of bed for splice sites", "<file>", NULL, &splicefilename);
   manopt(&optset, FLAG, 0, 'e', "edges", 
      "extract splice edges", NULL, NULL, &spliceedges);
  manopt(&optset, REQINTOPT, 0, 'V', "trans", 
      "generate random trans splicing events", NULL, NULL, &trans);
  manopt(&optset, FLAG, 0, 'D', "dist", 
      "generate dist splice sites", NULL, NULL, &distsplice);
  manopt(&optset, FLAG, 0, 'G', "simulate", 
      "simulate sequencing", NULL, NULL, &simulate);
  manopt(&optset, FLAG, 0, 'F', "isoforms", 
      "extract isoform sequences", NULL, NULL, &isoforms);
  manopt(&optset, REQUINTOPT, 0, 'C', "coverage", 
      "coverage for simulation", "<n>", NULL, &cov);


  unflagged = manopt_getopts(&optset, argc, argv);  
  if(unflagged->noofvalues > 1) {
    manopt_help(&optset, "unknown argument(s)\n");
  }


  assert(Pmis+Pins < 1.0);
  Pdel = 1.0 - Pmis - Pins;

  if(dbfilename) {
    fasta = bl_fastxGetSet(space, &dbfilename, 1, 1, 0, 0, 1);
  }

  if(splicefilename) {
    mytrack= bl_BEDread (space, splicefilename);
    geneset = bl_getGeneModelFromBEDtrack(space, mytrack);
  }

  alphabetlen = strlen(alphabet);
  if(readsfilename) { 
    readsdev = fopen(readsfilename, "w");
    if(readsdev == NULL) {
      fprintf(stderr, "couldn't open %s - exit forced", readsfilename);
      exit(-1);
    }
  }


  if(isoforms) {
    if(!dbfilename || !splicefilename)
      manopt_help(&optset, "database and splice annotation needed!\n");
    for(i=0; i < geneset->noofgenes; i++) { 
      sequence = bl_getGeneSequence(space, fasta, &geneset->genes[i]);
      fprintf(readsdev, ">%s\n%s\n", geneset->genes[i].id, sequence);
    }
  }

  if(simulate) { 
    if(!dbfilename || !splicefilename)
      manopt_help(&optset, "database and splice annotation needed!\n");
    bl_simulateGeneSetSequencing (space, readsdev, fasta, geneset, minreadlen, cov,
        alphabet, alphabetlen, minqual, maxqual, acc, Pmis, Pins);
  }

  if(trans) { 
    if(!splicefilename)
      manopt_help(&optset, "splice annotation needed!\n");
    if(distsplice) transtype = 'D';
    transset = bl_simulateTransSplicing (space, geneset, transtype, trans);
    transtrack = bl_getTrackFromGeneModel (space, transset);
    bl_BEDwrite(transtrack, dev); 
  }

  if(spliceedges) { 
    if(!splicefilename)
      manopt_help(&optset, "splice annotation needed!\n");
    bl_printSplicingEdges (space, dev, geneset);
  }

  if(subjectfilename) {  
    subjectdev = fopen(subjectfilename, "w");
    sequence = ALLOCMEMORY(space, NULL, char, reflen+1);
    fprintf(subjectdev, ">subject sequence (len: %d)\n", reflen);
    for(i=0; i < reflen; i++) {
      sequence[i]=alphabet[(Uint) RANDINT(alphabetlen)];
      fprintf(subjectdev,"%c", sequence[i]);
      if (i > 0 && i % 80 == 0) fprintf(subjectdev,"\n");
    }
    sequence[reflen] = '\0';
    fclose(subjectdev);
  }

  if(fiveprime) {
    fiveprimelen = strlen(fiveprime);
    fprintf(stderr,"5prime adapter: %s (%d)", fiveprime, fiveprimelen);
  }
  if(threeprime) {
    threeprimelen = strlen(threeprime);
    fprintf(stderr,"3prime adapter: %s (%d)", threeprime, threeprimelen);
  }


  if(subjectfilename && !matefilename && !split && !biseq) {
   
    bl_fastxPrintRandomReads(readsdev, sequence, reflen, 
      n, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
			     fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, NULL);
  }


  if(subjectfilename && splice && !matefilename && !biseq) {
    bl_fastxSimulateSpliceSites (space, sequence, reflen, 100, maxchildren, 
        alphabet, alphabetlen, acc, 
        Pmis, Pins, Pdel,
        minqual, maxqual, Pcis, mincisdist, maxcisdist, Pstrandswitch, 100);
  
  }

  if(subjectfilename && split && !matefilename && !biseq) {
    fprintf(stderr, "reflen: %d; maxreadlen %d\n", reflen, maxreadlen);
    bl_fastxPrintRandomSplitReads(readsdev, sequence, reflen, 
        n, minreadlen, maxreadlen, alphabet, alphabetlen, acc,
        Pmis, Pins, Pdel, fastq, minqual, maxqual,
        fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen);
  }

  if(subjectfilename && matefilename && !biseq) {
    matedev = fopen(matefilename, "w");

    bl_fastxPrintRandomMatePairs(readsdev, matedev, sequence, reflen, 
        n, minreadlen, maxreadlen, mindist, maxdist, 
        alphabet, alphabetlen, acc,
        Pmis, Pins, Pdel, fastq, minqual, maxqual,
        fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen);
    
    fclose(matedev);
  }

  if(subjectfilename && biseq){
    //assert(rate > 0);
    /* reverse complement of subject */
    revcomp = charDNAcomplement(space, sequence, reflen);
    
    /* convert nucleotides of each strand independently */
    for (i = 0; i < reflen; i++){
      if (sequence[i] == 'C' && RANDUNIT < rate) sequence[i] = 'T';
      if (revcomp[i] == 'C' && RANDUNIT < rate) revcomp[i] = 'T';      
    }

    /* Lister et al. protocol */
    if (biseq == 1){
      /* +FW reads */
      bl_fastxPrintRandomReads(readsdev, sequence, reflen, 
      (int)n/2, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "+FW");
      /* -FW reads */
      bl_fastxPrintRandomReads(readsdev, revcomp, reflen, 
      n-(int)n/2, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "-FW");
    }
    /* Cokus et al. protocol (+FW, +RC, -FW, -RC) */
    else {      
      /* +FW reads */
      bl_fastxPrintRandomReads(readsdev, sequence, reflen, 
      (int)n/4, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "+FW");
      /* -FW reads */
      bl_fastxPrintRandomReads(readsdev, revcomp, reflen, 
      (int)n/4, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "-FW");
      /* +RC reads */
      char *tmp = charDNAcomplement(space, sequence, reflen);      
      bl_fastxPrintRandomReads(readsdev, tmp, reflen, 
      (int)n/4, minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "+RC");
      FREEMEMORY(space, tmp);
      /* -RC reads */
      tmp = charDNAcomplement(space, revcomp, reflen);
      bl_fastxPrintRandomReads(readsdev, tmp, reflen, 
      n-3*((int)n/4), minreadlen, maxreadlen, alphabet, alphabetlen, acc, 
      Pmis, Pins, Pdel, fastq, minqual, maxqual, 
      fiveprime, fiveprimelen, threeprime, threeprimelen, polyAlen, "-RC");
      FREEMEMORY(space, tmp);
      FREEMEMORY(space, revcomp);
    }
  }
 
  if(readsfilename) fclose(readsdev);
  FREEMEMORY(space, sequence);
  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(space, unflagged);
  return 0;
}

