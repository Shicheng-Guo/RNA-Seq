#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file (can be stdin)"),
	make_option(c("-r", "--organism"), default="human", help="organsim name (default: %default)"),
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

#listMarts(host="www.ensembl.org")
#mart <- useMart("ensembl") -- obsolete
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
## use hg38 assembly
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")

#listDatasets(mart)
if(opt$organism=="human") {
    mart <- useDataset("hsapiens_gene_ensembl", mart)
} else {
    mart <- useDataset("mmusculus_gene_ensembl", mart)
}

#listAttributes(mart)
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
