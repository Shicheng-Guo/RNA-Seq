#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--countFiles"), help="input files containing read counts from HTSeq seperate by comma"),
  make_option(c("-f", "--gffFile"), help="input GFF file used as input to HTSeq"),
	make_option(c("-o", "--outDir"), default=".", help="output directory to keep results (default: %default)"),
  make_option(c("-p", "--pvalue"), default="0.05", help="adjusted p-value (default: %default)"),
  make_option(c("-t", "--treatment"), help="A comma seperated list about description of the data, example: WT,WT,KO,KO,WT,WT,KO,KO (must be equal to the number of input files)"),
	make_option(c("-c", "--condition"), help="A comma seperated list about additonal description of the data, example: CFUE,CFUE,CFUE,CFUE,PreMegE,PreMegE,PreMegE,PreMegE (must be equal to the number of input files)"),
  make_option(c("-r", "--processor"), default=1, help="number of processors to use (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$countFiles) | is.null(opt$gffFile) | is.null(opt$treatment)) {
	cat("\nProgram: dexseq.R (R script to identify differential expressed exons using DEXSeq)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(BiocParallel))

## start analysis
## next five lines to debug
#setwd("/Users/sachinpundhir/Lobster/project/rna-seq-analysis/matilda/result_CFUE/07_dexseq")
#opt <- NULL
#opt$countFiles <- "./Index2_tophat_out_Abundance.count,./Index4_tophat_out_Abundance.count,./Index5_tophat_out_Abundance.count,./Index6_tophat_out_Abundance.count,./Index7_tophat_out_Abundance.count,./Index12_tophat_out_Abundance.count"
#opt$gffFile <- "../04_cuffmerge/merged_known.gff"
#opt$treatment <- "WT,WT,WT,KO,KO,KO"
#opt$processor <- 50

## start analysis
countFiles <- unlist(strsplit(opt$countFiles, ","))
treatment <- unlist(strsplit(opt$treatment, ","))
sampleTable <- data.frame(row.names = countFiles, condition=treatment)

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= ~ sample + exon + condition:exon, flattenedfile=opt$gffFile)
BPPARAM = MulticoreParam(workers=opt$processor)
dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)
dxd = testForDEU(dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr = DEXSeqResults(dxd)

setwd(opt$outDir)
# plot quality measures
pdf("quality.pdf")
par(mfrow=c(2,2))
plotDispEsts(dxd)
plotMA(dxr, cex=0.8)
dev.off()

# print output as HTML file
DEXSeqHTML(dxr, FDR=opt$pvalue, color=c("#FF000080", "#0000FF80"))

## print differentially expressed exons
res <- as.data.frame(dxr)
res <- res[which(res$padj < opt$pvalue),c(1:21)]
write.table(res, "dexseq.xls", sep="\t", quote=F, row.names=F)

## save the session for further use
save.session("dexseq.session")
