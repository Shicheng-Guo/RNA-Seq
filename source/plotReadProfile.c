//#############################################################################################
//# A C program to plot read profile.                                                         #
//#############################################################################################

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<getopt.h>

#define PROGRAM "plotReadProfile (visualize read profile as block group, block and read arrangements)"
#define AUTHOR "University of Copenhagen, Denmark"
#define VERSION "1.0"
#define CONTACT "sachin@rth.dk"

char *quefile = NULL;
char *outDir = ".";
char *annfile = NULL;
char *ncrna   = NULL;
char *bgId    = NULL;
char *tmpName = "tmp";
char *progDir = "/home/users/sachin/software/myScripts/";

// temporary files
char *tmpNPRFile = NULL;
char *tmpBGProfileFile = NULL;
char *tmpBProfileFile = NULL;
char *tmpReadProfileFile = NULL;
char *tmpFinalProfileFile = NULL;

//Round Normalized Read count difference to two decimal places.
#define round(x) ((x)>0?(long)((x)+0.5):(long)((x)-0.5))

//===============
//  VARIABLES  // 
//===============

FILE *QUEFILE; //File handle for input query read profile(s) file.
FILE *NPRFILE; //File handle to output normalized profile file.
FILE *BG_PROFILE_FILE; //File handle to output read profile as block group arrangement
FILE *B_PROFILE_FILE; //File handle to out block read profile as line graph
FILE *READ_PROFILE_FILE; //File handle for output read profile as read arrangement
FILE *FINAL_FILE; //File handle for final output file
long int max_sread_x_profile, max_tread_x_profile; //Variable to store maximum start and total reads in Profile (X).
int maxXBox=0, maxYBox=0; //Variable to store the x and y end coordinate for the last plotted block.
int startX=1000, endX=0; //Variable to store the start and end coordinate for the whole block group while box and read plotting.
int genomicStart=0;

//================
//  STRUCTURES  //
//================

typedef struct {
	char *id;
	int start;
	int end;
	double height;
} tag_t;

//---------------------//
typedef struct {
	char *chrom;
	int start;
	int end;
	char *strand;
	double height;
	tag_t *tags;
	int nooftags;
} block_t;

//---------------------//
typedef struct {
	block_t *blocks;
	int noofblocks;
	char *id;
	char *chrom;
	int start;
	int end;
	char *strand;
	double expression;
	int tags;
	char *ncRNAid;
	char *ncRNAtype;
	char *ncRNAclass;
} cluster_t;

//---------------------//
typedef struct {
	cluster_t *clusters;
	int noofclusters;
} set_t;

//--------------------------------------------//

//===============
//  FUNCTIONS  //
//===============

void compProfile(int, cluster_t*); //Function to compute the normalized read profile.
void normalize(double*, double*, double*, long int, double, double); //Function to normalize the values of a read profile between 0 and 1.
void replace_special_char(char*); //Function to replace special characters in a string with '-';
void plotBox(int, cluster_t*); //Function to plot block groups in box format.

//--------------------------------------------//

void usage(void) {
	fprintf(stdout, "\nProgram: %s\n", PROGRAM);
	fprintf(stdout, "Author: %s\n", AUTHOR);
	fprintf(stdout, "Version %s\n", VERSION);
	fprintf(stdout, "Contact %s\n", CONTACT);
	fprintf(stdout, "\nUsage: \tplotReadProfile [OPTIONS]\n");
	fprintf(stdout, "\t-q <file>  [input file with mapped read profile(s)]\n");
	fprintf(stdout, "[OPTIONS]\n");
	fprintf(stdout, "\t-p <dir>   [program directory (default: /home/users/sachin/software/myScripts/)]\n");
	fprintf(stdout, "\t-a <file>  [annotation file with actual coordinates (optional)]\n");
	fprintf(stdout, "\t-o <dir>   [output directory to place plots (default: current)]\n");
	fprintf(stdout, "\t-n <str>   [plot for specific ncRNA id (optional)]\n\n");
	exit(1);
}

//--------------------------------------------//

void parseopt(int argc, char *argv[]) {
	char c;
	while((c = getopt(argc, argv, "q:p:a:o:n:b:")) != -1) {
		switch(c) {
			case 'q':
				quefile = optarg;
				break;
			case 'p':
				progDir = optarg;
				break;
			case 'a':
				annfile = optarg;
				break;
			case 'o':
				outDir = optarg;
				break;
			case 'n':
				ncrna = optarg;
				break;
			case 'b':
				bgId = optarg;
				break;
			case '?':
				//fprintf(stderr, "Invalid argument(s)...");
				usage();
				break;
		}
	}
}
//--------------------------------------------------------------------------------------------------//
//Function to read the annotation file and retrieve actual coordinates.
void retActCoor(char *filename, cluster_t *cluster1) {
	FILE *f; int found=0;
	//Open the file having annotations.
	if((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s\n", filename);
		exit(0);
	}

	int e; char c;
	while((c=getc(f)) != EOF) {
		char chrom[20], strand[10], id[50], count[10], source[50], class[50], type[50];
		int start, end;
		e = fscanf(f, "%s %d %d %s %s %s %s %s %s", chrom, &start, &end, id, count, strand, source, type, class);
		if((start >= cluster1->start && start <= cluster1->end) ||
		   (end >= cluster1->start && end <= cluster1->end) ||
		   (start >= cluster1->start && end <= cluster1->end) ||
		   (start <= cluster1->start && end >= cluster1->end)) {
			startX = start;
			endX = end;
			found = 1;
			break;
		}
	}

	if(found==0) {
		fprintf(stderr, "\nActual coordinates for the query are not present in %s\n", filename);
		//exit(0);
	}
	//printf("%d\t%d\n", startX, endX);
	fclose(f);
}

//--------------------------------------------------------------------------------------------------//
//Function to read the query and subject file with read profile(s).
set_t *read_cluster_file(char *filename) {
	set_t *set;
	set = calloc(1, sizeof(set_t));
	set->clusters = NULL;
	set->noofclusters = 0;

	FILE *f;
	//Open the file having read profile(s).
	if((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s\n", filename);
		exit(0);
	}

	int e=0; char c;
	while((c=getc(f)) != EOF) {
	//while(e != EOF) {
		char cluster_id[50], cluster_chrom[50], cluster_strand[50], tag_chrom[50], tag_strand[50], tag_id[50], dum[1000], ncRNA_id[50] = "", \
		ncRNA_type[50] = "", ncRNA_class[100] = "";
		int cluster_start, cluster_end, cluster_tagCount, cluster_blockCount, tag_start, tag_end, block_no;
		double cluster_expression, tag_height;

		//Check for header.
		int header = 0;
		//char c = getc(f);
		if(c == '>')
			header = 1;
		ungetc(c, f);

		if(header == 1) {
			//Parse information.
			e = fscanf(f, "%s %s %d %d %s %lf %d %d %s %s %s", cluster_id, cluster_chrom, &cluster_start, &cluster_end, cluster_strand, \
				&cluster_expression, &cluster_tagCount, &cluster_blockCount, ncRNA_id, ncRNA_type, ncRNA_class);

			//Create new cluster
			set->noofclusters++;
			set->clusters = realloc(set->clusters, (set->noofclusters+1)*sizeof(cluster_t));
			set->clusters[set->noofclusters].blocks = NULL;
			set->clusters[set->noofclusters].noofblocks = 0;

			//Assign Id.
			set->clusters[set->noofclusters].id = malloc(strlen(cluster_id)+1);
			strcpy(set->clusters[set->noofclusters].id, cluster_id);

			//Assign chromosome.
			set->clusters[set->noofclusters].chrom = malloc(strlen(cluster_chrom)+1);
			strcpy(set->clusters[set->noofclusters].chrom, cluster_chrom);

			//Assign start.
			set->clusters[set->noofclusters].start = cluster_start;

			//Assign end.
			set->clusters[set->noofclusters].end = cluster_end;

			//Assign strand.
			set->clusters[set->noofclusters].strand = malloc(strlen(cluster_strand)+1);
			strcpy(set->clusters[set->noofclusters].strand, cluster_strand);

			//Assign expression.
			set->clusters[set->noofclusters].expression = cluster_expression;

			//Assign ncRNA_id.
			set->clusters[set->noofclusters].ncRNAid = malloc(strlen(ncRNA_id)+1);
			strcpy(set->clusters[set->noofclusters].ncRNAid, ncRNA_id);

			//Assign ncRNAtype.
			set->clusters[set->noofclusters].ncRNAtype = malloc(strlen(ncRNA_type)+1);
			strcpy(set->clusters[set->noofclusters].ncRNAtype, ncRNA_type);

			//Assign ncRNAclass.
			set->clusters[set->noofclusters].ncRNAclass = malloc(strlen(ncRNA_class)+1);
			strcpy(set->clusters[set->noofclusters].ncRNAclass, ncRNA_class);
		}
		else {
			//Parse information.
			e = fscanf(f, "%s\t%d\t%d\t%s\t%lf\t%s\t%d\n", tag_chrom, &tag_start, &tag_end, tag_id, &tag_height, tag_strand, &block_no);

			//Get relative position of the tag.
			int newStart, newEnd;
			if(strcmp(tag_strand, "+")==0) { newStart = tag_start - cluster_start; newEnd = tag_end - cluster_start; }
			else { newStart = cluster_end - tag_end; newEnd = cluster_end - tag_start; }

			// get memory for new block, if new block starts
			if(block_no > set->clusters[set->noofclusters].noofblocks)
			{
				set->clusters[set->noofclusters].noofblocks++;
				set->clusters[set->noofclusters].blocks = realloc(set->clusters[set->noofclusters].blocks, \
					(set->clusters[set->noofclusters].noofblocks+1)*sizeof(block_t));
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags = 0;
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].start = newStart;
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].end = newEnd;
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].tags = NULL;
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].height = 0;

				// write block locus
				// chrom
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].chrom = malloc(strlen(tag_chrom)+1);
				strcpy(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].chrom, tag_chrom);
				// strand
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].strand = malloc(strlen(tag_strand)+1);
				strcpy(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].strand, tag_strand);
			}

			//Get memory for new tag.
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags++;
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].tags =
			realloc(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].tags,
				(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags+1)*sizeof(tag_t));

			//Assign tag locus Id.
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].
				tags[set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags].id = malloc(strlen(tag_id)+1);
			strcpy(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].
				tags[set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags].id, tag_id);

			//Assign tag start.
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].
				tags[set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags].start = newStart;
			
			//Assign tag end.
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].
				tags[set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags].end = newEnd;

			//Assign tag height.
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].
				tags[set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].nooftags].height = tag_height;

			//Assign block height
			set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].height += tag_height;

			//Change block boundaries 
			if(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].start > newStart)
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].start = newStart;
			if(set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].end < newEnd)
				set->clusters[set->noofclusters].blocks[set->clusters[set->noofclusters].noofblocks].end = newEnd;
		}
	}
	fclose(f);
	return set;
}
//--------------------------------------------------------------------------------------------------//
void compProfile(int a, cluster_t *cluster1) {
	int c, t, i, j;
	double sread_count, tread_count;
	double max_sread_x_profile=0, max_tread_x_profile=0, max_sread_y_profile=0, max_tread_y_profile=0;
	long int x_profile_len = 2*(cluster1->blocks[a].end - (cluster1->blocks[a].start-1));
	double x_profile[x_profile_len]; //Store Profile (X) in an array.
	double norm_x_profile[x_profile_len/2]; //Store Normalized Profile (X) in an array.
	double norm_sread_x_profile[x_profile_len/2]; //Store Normalized Start read count of Profile (X) in an array.
	double norm_tread_x_profile[x_profile_len/2]; //Store Normalized Total read count of Profile (X) in an array.

	//Compute read profile of cluster1.
	for(c=cluster1->blocks[a].start, i=0, j=0; c<=(cluster1->blocks[a].end); c++, j++) {
		sread_count=0; tread_count=0;
		for(t=1; t<=cluster1->blocks[a].nooftags; t++) {
			if(c>=cluster1->blocks[a].tags[t].start && c<=cluster1->blocks[a].tags[t].end) {
				tread_count += cluster1->blocks[a].tags[t].height;
			}
			if(c==cluster1->blocks[a].tags[t].start) {
				sread_count += cluster1->blocks[a].tags[t].height;
			}
		}
		//Put start read count information.
		x_profile[i] = sread_count;
		norm_sread_x_profile[j] = sread_count;
		if(sread_count > max_sread_x_profile) { max_sread_x_profile = sread_count; }
		i++;

		//Put total read count information.
		x_profile[i] = tread_count;
		norm_tread_x_profile[j] = tread_count;
		if(tread_count > max_tread_x_profile) { max_tread_x_profile = tread_count; }
		i++;
	}

	//Open output file for writing normalized profile aligment.
	if((NPRFILE = fopen(tmpNPRFile, "wt")) == NULL ) {
		fprintf(stderr, "Cannot open normalized profile alignment file for writing.\n");
	}

	//Normalize read profile.
	//normalize(norm_x_profile, norm_sread_x_profile, norm_tread_x_profile, x_profile_len/2, max_sread_x_profile, max_tread_x_profile);
	//normalize(norm_x_profile, norm_sread_x_profile, norm_tread_x_profile, x_profile_len/2, cluster1->blocks[a].height, cluster1->blocks[a].height);
	normalize(norm_x_profile, norm_sread_x_profile, norm_tread_x_profile, x_profile_len/2, cluster1->expression, cluster1->expression);

	fclose(NPRFILE);
}
//--------------------------------------------------------------------------------------------------//
//Function to normalize the values of a read profile between 0 and 1.
//void normalize(double *norm_profile, double *norm_sread_profile, double *norm_tread_profile, long int profile_len, int max_sread_profile, int max_tread_profile) {
void normalize(double *norm_profile, double *norm_sread_profile, double *norm_tread_profile, long int profile_len, double max_sread_profile, double max_tread_profile) {
	int i=0;
	while(i < profile_len) {
		//printf("%6.2f\t%6.2f\t%6.2f\t%6.2f\t", *(norm_sread_profile+i), max_sread_profile, *(norm_tread_profile+i), max_tread_profile);
		*(norm_sread_profile+i) = (double) round((*(norm_sread_profile+i)/max_sread_profile) * 100); //Normalized start read count.
		*(norm_tread_profile+i) = (double) round((*(norm_tread_profile+i)/max_tread_profile) * 100); //Normalized total read count.
		//printf("%6.2f\t%6.2f\t", *(norm_sread_profile+i), *(norm_tread_profile+i));
		*(norm_profile+i) = fabs(*(norm_sread_profile+i) - *(norm_tread_profile+i))/100;
		//printf("%6.2f\n", *(norm_profile+i));
		fprintf(NPRFILE, "%d\t%6.2f\n", i+1, *(norm_profile+i));
		i++;
	}
	//printf("\n");
}
//--------------------------------------------------------------------------------------------------//
//Function to replace special characters in a string with '-'.
void replace_special_char(char* string) {
	int i;
	for(i=0; i<strlen(string); i++) {
		if(string[i]==':' || string[i]=='\\' || string[i]=='/' || string[i]=='>' || \
		string[i]=='<' || string[i]=='\"' || string[i]=='\'' || string[i]=='|' || \
		string[i]=='\?' || string[i]=='*' || string[i]=='.' || string[i]=='(' || \
		string[i]==')' || string[i]=='_') {
			if(i==0) string[i]=' ';
			else string[i]='-';
		}
	}
}
//--------------------------------------------------------------------------------------------------//
void plotBox(int a, cluster_t *cluster1) {
	int mod_x1, mod_x2, mod_y1, mod_y2;
	char *id, *ncRNAid, *ncRNAtype, *ncRNAclass;
	int i; char shell_cmd[200];

	if(a==1) {
		id=malloc(100), ncRNAid=malloc(100), ncRNAtype=malloc(100), ncRNAclass=malloc(100);
		sprintf(id, "%s", cluster1->id); replace_special_char(id);
		sprintf(ncRNAid, "%s", cluster1->ncRNAid); replace_special_char(ncRNAid);
		sprintf(ncRNAtype, "%s", cluster1->ncRNAtype); replace_special_char(ncRNAtype);
		sprintf(ncRNAclass, "%s", cluster1->ncRNAclass); replace_special_char(ncRNAclass);
		fprintf(BG_PROFILE_FILE, "\\documentclass{standalone}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage{graphicx}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage{color}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage{tikz}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage{subfigure}\n");
		fprintf(BG_PROFILE_FILE, "\\usepackage{adjustbox}\n");
		fprintf(BG_PROFILE_FILE, "\\definecolor{dark-red}{RGB}{100,0,0}\n");
		fprintf(BG_PROFILE_FILE, "\\definecolor{dark-green}{RGB}{0,100,0}\n");
		fprintf(BG_PROFILE_FILE, "\\definecolor{dark-blue}{RGB}{0,0,100}\n");
		fprintf(BG_PROFILE_FILE, "\\definecolor{dark-yellow}{RGB}{255,185,15}\n");
		fprintf(BG_PROFILE_FILE, "\\begin{document}\n");
		fprintf(BG_PROFILE_FILE, "\\centering\n");
		fprintf(BG_PROFILE_FILE, "\\begin{adjustbox}{width=\\textwidth,height=\\textheight,keepaspectratio}\n");
		fprintf(BG_PROFILE_FILE, "\\begin{tikzpicture}[scale=1]\n");
		free(id);
		free(ncRNAid);
		free(ncRNAtype);
		free(ncRNAclass);
	}

	if(cluster1->blocks[a].start < maxXBox) {
		maxYBox=maxYBox+5;
		fprintf(BG_PROFILE_FILE, "\\draw [thick,white,fill=dark-red] (%d mm, %d mm) rectangle (%d mm, %d mm);", cluster1->blocks[a].start, maxYBox+5, cluster1->blocks[a].end, maxYBox+10);
		mod_x1 = cluster1->blocks[a].start; mod_x2 = (cluster1->blocks[a].start+cluster1->blocks[a].end)/2;
    mod_y1 = (maxYBox+5+(maxYBox+10))/2; mod_y2 = (maxYBox+5+(maxYBox+10))/2;
		fprintf(BG_PROFILE_FILE, "\\draw [white] (%d mm, %d mm)(%d mm, %d mm) node {Block %d};\n", mod_x1, mod_y1, mod_x2, mod_y2, a);
		if(cluster1->blocks[a].end > maxXBox) {
			maxXBox=cluster1->blocks[a].end;
		}
	}
	else {
		fprintf(BG_PROFILE_FILE, "\\draw [thick,white,fill=dark-red] (%d mm, 5 mm) rectangle (%d mm, 10 mm);", cluster1->blocks[a].start, cluster1->blocks[a].end);
		mod_x1 = cluster1->blocks[a].start; mod_x2 = (cluster1->blocks[a].start+cluster1->blocks[a].end)/2;
		mod_y1 = (5+10)/2; mod_y2 = (5+10)/2;
		fprintf(BG_PROFILE_FILE, "\\draw [white] (%d mm, %d mm)(%d mm, %d mm) node {Block %d};\n", mod_x1, mod_y1, mod_x2, mod_y2, a);
		if(cluster1->blocks[a].end > maxXBox) {
			maxXBox=cluster1->blocks[a].end;
		}
	}

	if(cluster1->blocks[a].start < startX) { startX = cluster1->blocks[a].start; }
	if(cluster1->blocks[a].end > endX) { endX = cluster1->blocks[a].end; }

	if(a==cluster1->noofblocks) {
		//fprintf(BG_PROFILE_FILE, "\\draw [thick, |-|] (%d mm, 3 mm) -- (%d mm, 3 mm);\n", startX, endX);
		fprintf(BG_PROFILE_FILE, "\\draw (%d mm, 0 mm) -- coordinate (x axis mid) (%d mm, 0 mm);\n", startX, endX);
		fprintf(BG_PROFILE_FILE, "\\foreach \\x in {%d,20,...,%d}\n", startX, endX);
		fprintf(BG_PROFILE_FILE, "\\draw (\\x mm,1 mm) -- (\\x mm,-2 mm)\n");
		fprintf(BG_PROFILE_FILE, "node[anchor=north] {};\n");

		//define genomic start position of block group
		if(strcmp(cluster1->strand, "+")==0) {	genomicStart = cluster1->start; }
		else { genomicStart = cluster1->end; }

		for(i=startX; i<=endX; i+=20) {
			if(i==startX) {
				fprintf(BG_PROFILE_FILE, "\\draw [black] (%d mm, -5 mm) (%d mm, -4 mm) node {%s:%d};\n", i-4, i-2, cluster1->chrom, genomicStart);
			}
			else {
				fprintf(BG_PROFILE_FILE, "\\draw [black] (%d mm, -5 mm) (%d mm, -4 mm) node {%d};\n", i, i, genomicStart);
			}
			if(strcmp(cluster1->strand, "+")==0) { genomicStart+=20; }
			else { genomicStart-=20; }
		}
		fprintf(BG_PROFILE_FILE, "\\draw [black] (-20 mm, 4 mm) (-15 mm, 4 mm) node {Block group};\n");
		fprintf(BG_PROFILE_FILE, "\\end{tikzpicture}\n");
		fprintf(BG_PROFILE_FILE, "\\end{adjustbox}\n");
		fprintf(BG_PROFILE_FILE, "\\end{document}\n");
		fclose(BG_PROFILE_FILE);

		sprintf(shell_cmd, "pdflatex %s > /dev/null", tmpBGProfileFile);
		system(shell_cmd);

	}
}
//--------------------------------------------------------------------------------------------------//
void plotReadProfile(cluster_t *cluster1) {
	int b, t, h;
	double y=5;
	char *id, *ncRNAid;
	int i; char shell_cmd[200];

	id=malloc(100), ncRNAid=malloc(100);
	sprintf(id, "%s", cluster1->id); id[0]='C';
	sprintf(ncRNAid, "%s", cluster1->ncRNAid); replace_special_char(ncRNAid);

	//Open output file for writing profile as line graph in tex file.
	if((READ_PROFILE_FILE = fopen(tmpReadProfileFile, "wt")) == NULL ) {
		fprintf(stderr, "Cannot open read profile file for writing.\n");
	}
	fprintf(READ_PROFILE_FILE, "\\documentclass{standalone}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage{graphicx}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage{color}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage{tikz}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage{subfigure}\n");
	fprintf(READ_PROFILE_FILE, "\\usepackage{adjustbox}\n");
	fprintf(READ_PROFILE_FILE, "\\definecolor{dark-red}{RGB}{100,0,0}\n");
	fprintf(READ_PROFILE_FILE, "\\definecolor{dark-green}{RGB}{0,100,0}\n");
	fprintf(READ_PROFILE_FILE, "\\definecolor{dark-blue}{RGB}{0,0,100}\n");
	fprintf(READ_PROFILE_FILE, "\\definecolor{dark-yellow}{RGB}{255,185,15}\n");
	fprintf(READ_PROFILE_FILE, "\\begin{document}\n");
	fprintf(READ_PROFILE_FILE, "\\centering\n");
	fprintf(READ_PROFILE_FILE, "\\begin{adjustbox}{width=\\textwidth,height=\\textheight,keepaspectratio}\n");
	fprintf(READ_PROFILE_FILE, "\\begin{tikzpicture}[scale=1]\n");

	//fprintf(READ_PROFILE_FILE, "\\normalsize \\textbf{Read\\hspace{1mm}profile:} (\\emph{Tag expression: \\color{dark-red}{$>=$10 reads}, \\color{dark-green}{$>=$1 and $<$10 reads}, \\color{dark-yellow}{$<$1 reads}})\n");
	
	//if(cluster1->expression <= 4000) {
	//Plot read profile of cluster1.
	for(b=1; b<=cluster1->noofblocks; b++) {
		for(t=1; t<=cluster1->blocks[b].nooftags; t++) {
			//for(h=0; h<cluster1->blocks[b].tags[t].height; h++) {
				if(cluster1->blocks[b].tags[t].height < 1)
					fprintf(READ_PROFILE_FILE, "\\draw [line width=0.05 mm, color=dark-yellow] (%d mm, %lf mm) -- (%d mm, %lf mm);\n", cluster1->blocks[b].tags[t].start, y, cluster1->blocks[b].tags[t].end, y);
				else if(cluster1->blocks[b].tags[t].height >= 1 && cluster1->blocks[b].tags[t].height < 5)
					fprintf(READ_PROFILE_FILE, "\\draw [line width=0.10 mm, color=dark-green] (%d mm, %lf mm) -- (%d mm, %lf mm);\n", cluster1->blocks[b].tags[t].start, y, cluster1->blocks[b].tags[t].end, y);
				else if(cluster1->blocks[b].tags[t].height >= 5 && cluster1->blocks[b].tags[t].height < 10)
					fprintf(READ_PROFILE_FILE, "\\draw [line width=0.15 mm, color=dark-green] (%d mm, %lf mm) -- (%d mm, %lf mm);\n", cluster1->blocks[b].tags[t].start, y, cluster1->blocks[b].tags[t].end, y);
				else if(cluster1->blocks[b].tags[t].height >= 10 && cluster1->blocks[b].tags[t].height < 20)
					fprintf(READ_PROFILE_FILE, "\\draw [line width=0.20 mm, color=dark-red] (%d mm, %lf mm) -- (%d mm, %lf mm);\n", cluster1->blocks[b].tags[t].start, y, cluster1->blocks[b].tags[t].end, y);
				else if(cluster1->blocks[b].tags[t].height >= 20)
					fprintf(READ_PROFILE_FILE, "\\draw [line width=0.25 mm, color=dark-red] (%d mm, %lf mm) -- (%d mm, %lf mm);\n", cluster1->blocks[b].tags[t].start, y, cluster1->blocks[b].tags[t].end, y);
			//}
			y = y + 0.25;
		}
		fprintf(READ_PROFILE_FILE, "\\draw [black] (%d mm, %lf mm)(%d mm, %lf mm) node {%0.2f};\n", cluster1->blocks[b].start, y+2, (cluster1->blocks[b].start+cluster1->blocks[b].end)/2, y+2, cluster1->blocks[b].height);
		//printf("\\draw [black] (%d mm, %lf mm)(%d mm, %lf mm) node {%lf};\n", cluster1->blocks[b].start, y, cluster1->blocks[b].end, y, cluster1->blocks[b].height);
		y = y + 5;
	}//}
	//else {
	//	fprintf(READ_PROFILE_FILE, "\\draw [black] (0 mm, 10 mm) (0 mm, 10 mm) node {{\\color{dark-blue}{Total reads in query block group are $>$4000. Too large to accommodate in one image.}}};\n");
	//}

	fprintf(READ_PROFILE_FILE, "\\draw (0 mm, 0 mm) -- coordinate (x axis mid) (%d mm, 0 mm);\n", cluster1->end-cluster1->start);
	fprintf(READ_PROFILE_FILE, "\\foreach \\x in {0,20,...,%d}\n", cluster1->end-cluster1->start);
	fprintf(READ_PROFILE_FILE, "\\draw (\\x mm,1 mm) -- (\\x mm,-2 mm)\n");
	fprintf(READ_PROFILE_FILE, "node[anchor=north] {};\n");

	//define genomic start position of block group
	if(strcmp(cluster1->strand, "+")==0) {	genomicStart = cluster1->start; }
	else { genomicStart = cluster1->end; }

	for(i=startX; i<=endX; i+=20) {
		if(i==startX) {
			fprintf(READ_PROFILE_FILE, "\\draw [black] (%d mm, -5 mm) (%d mm, -4 mm) node {%s:%d};\n", i-4, i-2, cluster1->chrom, genomicStart);
		}
		else {
			fprintf(READ_PROFILE_FILE, "\\draw [black] (%d mm, -5 mm) (%d mm, -4 mm) node {%d};\n", i, i, genomicStart);
		}
		if(strcmp(cluster1->strand, "+")==0) { genomicStart+=20; }
		else { genomicStart-=20; }
	}
	fprintf(BG_PROFILE_FILE, "\\draw [black] (-20 mm, 4 mm) (-15 mm, 4 mm) node {Read profile};\n");
	//fprintf(READ_PROFILE_FILE, "\\draw [thick, |-|, color=dark-red] (0 mm, -6 mm) -- (%d mm, -6 mm);\n", cluster1->end-cluster1->start);
	fprintf(READ_PROFILE_FILE, "\\end{tikzpicture}\n");
	fprintf(READ_PROFILE_FILE, "\\end{adjustbox}\n");
	fprintf(READ_PROFILE_FILE, "\\end{document}\n");
	fclose(READ_PROFILE_FILE);

	sprintf(shell_cmd, "pdflatex %s > /dev/null", tmpReadProfileFile);
	system(shell_cmd);

	free(id);
	free(ncRNAid);
}
//--------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
	//Parse and validate the arguments passed to DeepAlign.
	parseopt(argc, argv);
	if(argc < 3) { usage(); }

	//Read Query read profile(s).
	set_t *querySet = read_cluster_file(quefile);

	//Compare all cluster(s) in query file against all cluster(s) in subject file.
	cluster_t *query;
	int q, qb;
	char shell_cmd[500];
	char *bgIdF, *query_id, *query_ncRNAid;
	int x_axis_len;
	bgIdF = malloc(200);
	query_id = malloc(200);
	query_ncRNAid = malloc(200);
	tmpNPRFile=malloc(200);
	tmpBGProfileFile=malloc(200);
	tmpBProfileFile=malloc(200);
	tmpReadProfileFile=malloc(200);
	tmpFinalProfileFile=malloc(200);
	tmpName=malloc(200);
	
	//Run through all blocks in the query read profile(s).
	for(q=1; q<=querySet->noofclusters; q++) {
		//Continue if ncRNA name does not match to input criteria
		if(ncrna!=NULL && strcmp(querySet->clusters[q].ncRNAid, ncrna)!=0) { continue; }
		sprintf(bgIdF, ">%s", bgId);
		if(bgId!=NULL && strcmp(querySet->clusters[q].id, bgIdF)!=0) { continue; }
		//if(querySet->clusters[q].noofblocks<2) { continue; }

		// defined name for the output files
		sprintf(tmpName, "%s_%d_%d", querySet->clusters[q].chrom, querySet->clusters[q].start, querySet->clusters[q].end);
		// initialze temporary file names
		sprintf(tmpNPRFile, "%s.npr", tmpName);
		sprintf(tmpBGProfileFile, "%s.bgp.tex", tmpName);
		sprintf(tmpBProfileFile, "%s.bp.tex", tmpName);
		sprintf(tmpReadProfileFile, "%s.rp.tex", tmpName);
		sprintf(tmpFinalProfileFile, "%s.tex", tmpName);

		//Open output file for writing profile as block group arrangement in tex file.
		if((BG_PROFILE_FILE = fopen(tmpBGProfileFile, "wt")) == NULL ) {
			fprintf(stderr, "Cannot open block alignment file for writing.\n");
		}
		maxXBox=0, maxYBox=0, startX=1000, endX=0;
		for(qb = 1; qb <= querySet->clusters[q].noofblocks; qb++) {
			//Set cluster q of the query cluster(s) as query.
			query = &querySet->clusters[q];
	
			x_axis_len = (query->blocks[qb].end - (query->blocks[qb].start -1));
			sprintf(query_id, "%s", query->id);
			query_id[0]='C';
			
			sprintf(query_ncRNAid, "%s", query->ncRNAid);

			//Substitute special characters in query_ncRNAid with '-'.
			replace_special_char(query_ncRNAid);

			//Only plot profiles for clusters <= threshold
			//if(x_axis_len > 40) { continue; }
		
			//Determine the normalized read profile pattern of query profile.
			compProfile(qb, query);

			//Arguments: 'PNG file', 'x-axis end', 'output dir'
			sprintf(shell_cmd, "R --no-save --vanilla --slave < %s/plotblockProfile.r --args %s %s_%s_%d.pdf %d %s_%s_%s_%s_%d", progDir, tmpNPRFile, query_id, query_ncRNAid, qb, x_axis_len, query_id, query_ncRNAid, query->ncRNAtype, query->ncRNAclass, qb);
			system((char*)shell_cmd);
			//printf("%s\n", shell_cmd);

			//Print program status in STDOUT.
			printf("%s_%s_%d\tDone\n", query->id, query->ncRNAid, qb);
		
			//Plot block group in box alignment format.
			plotBox(qb, query);
		}

		//Retrieve actual coordinates
		if(annfile!=NULL) {	retActCoor(annfile, query); }

		//Plot read profile.
		plotReadProfile(query);

		//Open output file for writing final profile in tex file.
		if((FINAL_FILE = fopen(tmpFinalProfileFile, "wt")) == NULL ) {
			fprintf(stderr, "Cannot open final profile file for writing.\n");
		}

		//create final output file for read profiles
		fprintf(FINAL_FILE, "\\documentclass[a4paper, 10pt]{article}\n");
		fprintf(FINAL_FILE, "\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}\n");
		fprintf(FINAL_FILE, "\\usepackage{graphicx}\n");
		fprintf(FINAL_FILE, "\\usepackage{color}\n");
		fprintf(FINAL_FILE, "\\usepackage{tikz}\n");
		fprintf(FINAL_FILE, "\\usepackage{subfigure}\n");
		fprintf(FINAL_FILE, "\\usepackage{adjustbox}\n");
		fprintf(FINAL_FILE, "\\definecolor{dark-red}{RGB}{100,0,0}\n");
		fprintf(FINAL_FILE, "\\definecolor{dark-green}{RGB}{0,100,0}\n");
		fprintf(FINAL_FILE, "\\definecolor{dark-blue}{RGB}{0,0,100}\n");
		fprintf(FINAL_FILE, "\\begin{document}\n");
		fprintf(FINAL_FILE, "\\begin{center}\n");
		fprintf(FINAL_FILE, "\\Large{\\textbf{\\underline{deepBlockAlign\\hspace{1mm}v1.0}}}\n");
		fprintf(FINAL_FILE, "\\end{center}\n");
		fprintf(FINAL_FILE, "\\vspace{0.3cm}\n");
		fprintf(FINAL_FILE, "\\begin{center}\n");
		fprintf(FINAL_FILE, "\\begin{adjustbox}{width=\\textwidth,height=\\textheight,keepaspectratio}\n");
		fprintf(FINAL_FILE, "\\includegraphics[type=pdf,ext=.pdf,read=.pdf,scale=1]{%s.bgp}\n", tmpName);
		fprintf(FINAL_FILE, "\\end{adjustbox}\n");
		fprintf(FINAL_FILE, "\\end{center}\n");
		fprintf(FINAL_FILE, "\\vspace{0.3cm}\n");
		fprintf(FINAL_FILE, "\\begin{center}\n");
		int i;
		for(i=1; i<=querySet->clusters[q].noofblocks; i++) { 
			fprintf(READ_PROFILE_FILE, "\\includegraphics[scale=1]{%s_%s_%d.pdf}\n", query_id, query_ncRNAid, i);
		}
		fprintf(FINAL_FILE, "\\end{center}\n");
		fprintf(FINAL_FILE, "\\vspace{0.3cm}\n");
		fprintf(FINAL_FILE, "\\begin{center}\n");
		fprintf(FINAL_FILE, "\\begin{adjustbox}{width=\\textwidth,height=\\textheight,keepaspectratio}\n");
		fprintf(FINAL_FILE, "\\includegraphics[type=pdf,ext=.pdf,read=.pdf,scale=1]{%s.rp}\n", tmpName);
		fprintf(FINAL_FILE, "\\end{adjustbox}\n");
		fprintf(FINAL_FILE, "\\end{center}\n");
		fprintf(FINAL_FILE, "\\end{document}\n");
		fclose(FINAL_FILE);

		sprintf(shell_cmd, "pdflatex %s > /dev/null", tmpFinalProfileFile);
		system(shell_cmd);

		// move all files to output directory
		if(outDir!=".") {
			sprintf(shell_cmd, "mv %s* %s", tmpName, outDir);
			system(shell_cmd);
			sprintf(shell_cmd, "mv %s_%s* %s", query_id, query_ncRNAid, outDir);
			system(shell_cmd);
		}

		//Delete intermediate files.
		printf("\nDelete intermediate files.....\n");
		printf("\nDone\n");
		sprintf(shell_cmd, "rm %s/*.log", outDir);
		system(shell_cmd);
		sprintf(shell_cmd, "rm %s/*.aux", outDir);
		system(shell_cmd);
		sprintf(shell_cmd, "rm %s/*.tex", outDir);
		system(shell_cmd);
		sprintf(shell_cmd, "rm %s/*.npr", outDir);
		system(shell_cmd);
		for(i=1; i<=querySet->clusters[q].noofblocks; i++) { 
			sprintf(shell_cmd, "rm %s/%s_%s_%d.pdf", outDir, query_id, query_ncRNAid, i);
			system(shell_cmd);
		}
	}
	free(query_id);
	free(query_ncRNAid);
	free(tmpNPRFile);
	free(tmpBGProfileFile);
	free(tmpBProfileFile);
	free(tmpReadProfileFile);
	free(tmpFinalProfileFile);

	return 0;
}
//--------------------------------------------------------------------------------------------------//
