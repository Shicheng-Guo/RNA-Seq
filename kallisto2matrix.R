#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--configFile"), help="configuration file containing kallisto mapped file information (format: <id> <dir> | can be stdin)"),
	make_option(c("-o", "--outFile"), help="output file name (name.gene_counts and name.transcript_counts will be created)"),
	make_option(c("-g", "--genome"), default="mm9", help="organism genome (default: %default; Options: mm9, mm10, hg19, hg38)"),
	make_option(c("-t", "--tpm"), help="create tpm matrix (default: count)", action="store_true"),
	make_option(c("-s", "--sleuth"), help="create count matrix after batch correction using sleuth (default: count)", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$configFile) | is.null(opt$outFile)) {
	cat("\nProgram: kallisto2matrix.R (R script to compute gene count matrix from kallisto results)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

## read configuration file
if(identical(opt$configFile, "stdin")==T) {
    configFile <- read.table(file("stdin"), header=F)
} else {
    configFile <- read.table(opt$configFile, header=F)
}

## NOTE: external_gene_name is only available with host=grch37.ensembl.org

if(opt$genome=="hg19") {
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm9") {
    #mart = useMart(host = "may2012.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else if(opt$genome=="hg38") {
    mart = useMart(host = "dec2017.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm10") {
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else {
    cat("Unknown genome provided\n");
    print_help(parser)
    q()
}

## read each file into a data frame
abundance.lst <- apply(configFile, 1, function(x) read.table(sprintf("%s/abundance.tsv", x[2]), header=T))
abundance.df <- do.call(cbind.data.frame, abundance.lst)

if(is.null(opt$tpm)) { 
    ## retrieve columns containing counts
    est_counts <- abundance.df[,c(which(names(abundance.df) %in% "est_counts"))]
    est_counts <- cbind(abundance.df$target_id, est_counts)
    colnames(est_counts) <- c("target_id", as.vector(configFile$V1))

    write.table(est_counts, sprintf("%s.transcript_counts", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
} else {
    ## retrieve columns containing tpm
    tpm <- abundance.df[,c(which(names(abundance.df) %in% "tpm"))]
    tpm <- cbind(abundance.df$target_id, tpm)
    colnames(tpm) <- c("target_id", as.vector(configFile$V1))

    write.table(tpm, sprintf("%s.transcript_tpm", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
}

if(is.null(opt$tpm)) { 
    ## determine external gene name for ensembl id
    tx2gene <- getBM(filters="ensembl_transcript_id", attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), values=est_counts$target_id, mart=mart)
    files <- apply(configFile, 1, function(x) sprintf("%s/abundance.tsv", x[2]))

    ## collapse transcript counts into gene counts
    gene_counts <- tximport(files, type="kallisto", tx2gene=tx2gene[,c(1,3)])
    gene_counts <- cbind(as.data.frame(rownames(gene_counts$counts)), gene_counts$counts)
    colnames(gene_counts) <- c("target_id", as.vector(configFile$V1))

    write.table(gene_counts, sprintf("%s.gene_counts", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
} else {
    ## determine external gene name for ensembl id
    tx2gene <- getBM(filters="ensembl_transcript_id", attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), values=tpm$target_id, mart=mart)
    files <- apply(configFile, 1, function(x) sprintf("%s/abundance.tsv", x[2]))

    ## collapse transcript tpm into gene tpm
    gene_tpm <- tximport(files, type="kallisto", tx2gene=tx2gene[,c(1,3)], countsFromAbundance="lengthScaledTPM")
    gene_tpm <- cbind(as.data.frame(rownames(gene_tpm$counts)), gene_tpm$counts)
    colnames(gene_tpm) <- c("target_id", as.vector(configFile$V1))

    write.table(gene_tpm, sprintf("%s.gene_tpm", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
}

q()
