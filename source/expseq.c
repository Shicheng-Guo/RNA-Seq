#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "libs/biofilesio.h"

/******************
 * FUNCTIONS
 *****************/

void writeHeader(char *reads_File, char *annotations_File, FILE *file)
{	
	time_t t;
	t = time(NULL);
	if(file != NULL){
		fprintf(file,"# expseq.x started at %s", ctime(&t));
		fprintf(file,"# mapping file: %s\n", reads_File);
		fprintf(file,"# annotation file: %s\n", annotations_File);
		fprintf(file,"");	
	}else{
		printf("# expseq.x started at %s", ctime(&t));
		printf("# mapping file: %s\n", reads_File);
		printf("# annotation file: %s\n", annotations_File);
		printf("");
	}
}

void printUsage(void)
{ 
	printf("\n");
	printf("==============================================\n");
	printf("                  ExpSeq                      \n");
	printf("analyze gene expression using RNAseq reads\n");
	printf("==============================================\n\n");
	printf("usage: expseq.x -r <file> -a <file> [ARGUMENTS]\n\n");
	printf("[INPUT]\n");
	printf(" -r  <file>    file with mapped reads (segemehl output)\n");
	printf(" -a  <file>    file with annotated genes (bed format)\n");
	printf(" -n  <int>     normalize read expression based on:\n");
	printf("                1) multiple hits (default)\n");
	printf("                2) multiple hits + mappable reads\n");
	printf("                3) multiple hits + all reads\n");
	printf("                4) off\n");
	printf(" -c  <double>  expression cutoff (default = 0)\n");
	printf(" -o  <file>    output file\n");
	printf("[ARGUMENTS]\n");
	printf(" -h           this helpfull message\n");
	printf("[VERSION]\n");
	printf(" 11-16-2010\n");	
}

geneset_t *getExpression(tagset_t *tagset, geneset_t *geneset, int normalize, double expressioncutoff, FILE *file)
/* GET EXPRESSION */
{
	int i;
	for(i = 1; i <= geneset->noofgenes; i++){
		int j; 
		if(geneset->genes[i].exonCount == 0){
			for(j = 1; j <= tagset->nooftags; j++){
				if(strcmp(tagset->tags[j].chrom, geneset->genes[i].chrom)!=0){continue;} 
				if(strcmp(tagset->tags[j].strand, geneset->genes[i].strand)!=0){continue;}
			
				if((tagset->tags[j].start > geneset->genes[i].start && tagset->tags[j].start < geneset->genes[i].end) || 
				   (tagset->tags[j].end > geneset->genes[i].start && tagset->tags[j].end < geneset->genes[i].end) ||
				   (tagset->tags[j].start < geneset->genes[i].start && tagset->tags[j].end > geneset->genes[i].end)){
					int x;
					for(x = tagset->tags[j].start; x <= tagset->tags[j].end; x++){
						if(x >= geneset->genes[i].start && x <= geneset->genes[i].end){
							if(normalize == 1)
								geneset->genes[i].expression += (tagset->tags[j].expression / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
							else if(normalize == 2)
								geneset->genes[i].expression += ((tagset->tags[j].expression / tagset->noofreads) / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
							else if(normalize == 3)
								geneset->genes[i].expression += ((tagset->tags[j].expression / tagset->noofallreads) / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
							else if(normalize == 4)
								geneset->genes[i].expression += ((tagset->tags[j].noofreads / (double) (tagset->tags[j].end - tagset->tags[j].start + 1)));
						}
					}
				}
			}
		}
		else{
			for(j = 1; j <= tagset->nooftags; j++){
				if(strcmp(tagset->tags[j].chrom, geneset->genes[i].chrom)!=0){continue;} 
				if(strcmp(tagset->tags[j].strand, geneset->genes[i].strand)!=0){continue;}
				int e;
				for(e = 0; e <= geneset->genes[i].noofexons; e++){
					if((tagset->tags[j].start > geneset->genes[i].exons[e].start && tagset->tags[j].start < geneset->genes[i].exons[e].end) || 
					   (tagset->tags[j].end > geneset->genes[i].exons[e].start && tagset->tags[j].end < geneset->genes[i].exons[e].end) ||
					   (tagset->tags[j].start < geneset->genes[i].exons[e].start && tagset->tags[j].end > geneset->genes[i].exons[e].end)){
						int x;
						for(x = tagset->tags[j].start; x <= tagset->tags[j].end; x++){
							if(x >= geneset->genes[i].exons[e].start && x <= geneset->genes[i].exons[e].end){
								if(normalize == 1)
									geneset->genes[i].expression += (tagset->tags[j].expression / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
								else if(normalize == 2)
									geneset->genes[i].expression += ((tagset->tags[j].expression / tagset->noofreads) / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
								else if(normalize == 3)
									geneset->genes[i].expression += ((tagset->tags[j].expression / tagset->noofallreads) / (double) (tagset->tags[j].end - tagset->tags[j].start + 1));
								else if(normalize == 4)
									geneset->genes[i].expression += ((tagset->tags[j].noofreads / (double) (tagset->tags[j].end - tagset->tags[j].start + 1)));
							}
						}
					}
				}
			}
		}
		if((geneset->genes[i].expression / ((double)geneset->genes[i].end - (double)geneset->genes[i].start + 1)) >= expressioncutoff){
			if(normalize == 1 || normalize == 4)
				fprintf(file,"%s\t%lf\n",geneset->genes[i].id,(geneset->genes[i].expression / ((double)geneset->genes[i].end - (double)geneset->genes[i].start + 1)));
			else 
				fprintf(file,"%s\t%E\n",geneset->genes[i].id,(geneset->genes[i].expression / ((double)geneset->genes[i].end - (double)geneset->genes[i].start + 1)));
		}
	}
	return geneset;
}


/* MAIN FUNCTION */
int main(int argc, char *argv[])
{  
	// check user input
	if((argc) == 1){printUsage(); exit(0);}
	char *reads_File = NULL;
	char *annotations_File = NULL;
	char *result_File = NULL;
	int normalize = 1;
	double expressioncutoff = 0;
	
	int i;
	for (i = 1; i < argc; i++)  
    {	
		if (strcmp(argv[i], "-r") == 0)
			reads_File = argv[i+1];
		
		if (strcmp(argv[i], "-a") == 0)
		    annotations_File = argv[i+1];
		
		if (strcmp(argv[i], "-n") == 0)
		    normalize = atoi(argv[i+1]);
		
		if (strcmp(argv[i], "-c") == 0)
		    expressioncutoff = atoi(argv[i+1]);
		
		if (strcmp(argv[i], "-o") == 0)
		    result_File = argv[i+1];
		
		if (strcmp(argv[i], "-h") == 0){
			printUsage(); exit(0);}
	}

	if(!reads_File || !annotations_File)
	{
		printUsage(); exit(0);
	}
	
	// open file for results 
	FILE *outFile = stdout;
	if(result_File != NULL)
		outFile = fopen(result_File,"w");
	
	// write header
	writeHeader(reads_File, annotations_File, outFile);	

	// read input files
	tagset_t *myTags = read_segemehl_file(reads_File);
	geneset_t *myAnnotations = read_bed_file(annotations_File);
	
	// get expression
	myAnnotations = getExpression(myTags, myAnnotations, normalize, expressioncutoff, outFile);
	
	fclose(outFile);
	
	return(0);
}
