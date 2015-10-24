#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file"),
	make_option(c("-r", "--organism"), default="human", help="organsim name (default: %default)"),
	make_option(c("-d", "--header"), default=T, help="if file has header (default: %default)"),
	make_option(c("-t", "--tab"), help="file is tab separated", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: convertId.R (R script to convert Ids using BioMart)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

## read input file
if(is.null(opt$tab)) {
    data <- read.table(opt$inFile, sep=",", header=as.logical(opt$header))
} else {
    data <- read.table(opt$inFile, header=as.logical(opt$header))
}

#listMarts()
mart <- useMart("ensembl")

#listDatasets(mart)
if(opt$organism=="human") {
    mart <- useDataset("hsapiens_gene_ensembl", mart)
} else {
    mart <- useDataset("mmusculus_gene_ensembl", mart)
}

#listAttributes(mart)
result <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name"), values=data[,1], mart=mart)
#data$geneName <- as.vector(unlist(apply(data, 1, function(x) result[which(result$ensembl_gene_id==x[1]),2])))

result <- merge(data, result, by.x=colnames(data)[1], by.y="ensembl_gene_id")

outfile=sprintf("%s.geneId", opt$inFile)
write.table(result, outfile, sep="\t", quote=F, row.names=F, col.names=as.logical(opt$header))
