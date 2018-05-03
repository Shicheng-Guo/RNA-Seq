#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file (can be stdin)"),
	make_option(c("-r", "--organism"), default="hg19", help="organsim name (hg19, mm9, hg38, mm10) (default: %default)"),
	make_option(c("-d", "--header"), default=T, help="if file has header (default: %default)"),
	make_option(c("-t", "--tab"), help="file is tab separated", action="store_true"),
	make_option(c("-b", "--bed"), help="also include gene coordinates in output", action="store_true")
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
    if(identical(opt$inFile, "stdin")==T) {
        data <- read.table(file("stdin"), sep=",", header=as.logical(opt$header))
    } else {
        data <- read.table(opt$inFile, sep=",", header=as.logical(opt$header))
    }
} else {
    if(identical(opt$inFile, "stdin")==T) {
        data <- read.table(file("stdin"), header=as.logical(opt$header))
    } else {
        data <- read.table(opt$inFile, header=as.logical(opt$header))
    }
}

data[,1] <- gsub("\\..*", "", data[,1])

################################################
## BIOMART USAGE
# listMarts(); ensembl = useMart("ENSEMBL_MART_ENSEMBL"); listDatasets(ensembl); listAttributes(mart)
# Use the http://www.ensembl.org website and go down the bottom of the page. Click on ’view in Archive’ and select the archive you need. Copy the url and use that url as shown below to connect to the specified BioMart database. The example below shows how to query Ensembl 54.
# > listMarts(host = "may2009.archive.ensembl.org")
# > ensembl54 = useMart(host = "may2009.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
################################################

## NOTE: external_gene_name is only available with host=grch37.ensembl.org

if(opt$organism=="hg19") {
    #mart = useMart(host = "feb2014.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$organism=="mm9") {
    mart = useMart(host = "may2012.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else if(opt$organism=="hg38") {
    mart = useMart(host = "dec2017.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$organism=="mm10") {
    #mart = useMart(host = "dec2017.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else {
	cat("Unknown organism provided\n");
	print_help(parser)
	q()
}

if(is.null(opt$bed)) {
    result <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name"), values=data[,1], mart=mart)
    #data$geneName <- as.vector(unlist(apply(data, 1, function(x) result[which(result$ensembl_gene_id==x[1]),2])))

    result <- merge(data, result, by.x=colnames(data)[1], by.y="ensembl_gene_id")

    outfile=sprintf("%s.geneId", opt$inFile)
    write.table(result, "", sep="\t", quote=F, row.names=F, col.names=as.logical(opt$header))
} else {
    result <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand"), values=data[,1], mart=mart)
    #data$geneName <- as.vector(unlist(apply(data, 1, function(x) result[which(result$ensembl_gene_id==x[1]),2])))

    result <- merge(data, result, by.x=colnames(data)[1], by.y="ensembl_gene_id")

    outfile=sprintf("%s.bed", opt$inFile)
    start <- ncol(result)-3
    result[,start] <- sprintf("chr%s", result[,start])
    end <- ncol(result)
    result <- result[,c(start:end,1:(start-1))]
    result[which(result$strand=="1"),]$strand <- "+"
    result[which(result$strand=="-1"),]$strand <- "-"
    write.table(result, "", sep="\t", quote=F, row.names=F, col.names=as.logical(opt$header))
}
