## EXAMPLE RUN
## R --no-save --vanilla --slave < seqbias.r --args 0 BloodGm12878Rep1.bam.yaml /home/projects/deepalign/data/genomes/hg19/hg19.fa BloodGm12878Rep1.bam
## R --no-save --vanilla --slave < seqbias.r --args 0 BloodGm12878Rep1.bam.yaml /home/projects/deepalign/data/genomes/hg19/hg19.fa DP.bed

## ARGUMENTS
## ARGV[1]: option
## ARGV[2]:
##          option1: output yaml model file
##          option2: input directory having .yaml files
## ARGV[3]: input fasta file (index)
## ARGV[4]: 
##          option 1: input bam file (index)
##          option 2: input bed file

library("seqbias")

## load arguments
ARGV <- commandArgs(TRUE)

## mode 1: create models; mode 2: predict bias using model
if(ARGV[1]==0) {
	yaml_file <- ARGV[2]
	fasta_file <- ARGV[3]
	bam_file <- ARGV[4]

	sb <- seqbias.fit(fasta_file, bam_file, L=15, R=15)
	seqbias.save(sb, yaml_file)
} else if(ARGV[1]==1) {
	yaml_dir <- ARGV[2]
	fasta_file <- ARGV[3]
	coor_file <- ARGV[4]
	#yaml_dir <- "."
	#fasta_file <- "/home/projects/deepalign/data/genomes/hg19/hg19.fa"
	#coor_file <- "test"
	ylist <- list.files(pattern = paste(yaml_dir, "*.yaml$", sep="/"))
	clist <- read.table(coor_file)

	## output euclidean distance for each coordinate
	sink(ARGV[5])

	for(i in 1:length(clist[,1])) {
		coor <- clist[i,]
		mat <- data.frame()
		for(j in 1:length(ylist)) {
			sb <- seqbias.load(fasta_file, ylist[j])
			I <- GRanges(coor$V1, IRanges( coor$V2, coor$V3), strand = coor$V4)
			bias <- seqbias.predict(sb, I)
			mat <- rbind(mat, bias$`1`)
		}
		bias <- as.vector(dist(mat, method="euclidean"))
		#bias <- as.vector(apply(mat, 1, sd))
		cat(bias)
		## uncomment, if need bias for each coordinate in new line
		#cat("\n")
	}
	sink()
} else if(ARGV[1]==2) {
	yaml_file <- ARGV[2]
	fasta_file <- ARGV[3]
	coor <- unlist(strsplit(ARGV[4], "_"))

	sb <- seqbias.load(fasta_file, yaml_file)
	I <- GRanges(coor[1], IRanges(as.numeric(coor[2]), as.numeric(coor[3])), strand = coor[4])
	bias <- seqbias.predict(sb, I)
	cat(bias$`1`)
}
