#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--matrixFile"), help="input read count file"),
    make_option(c("-t", "--treatment"), help="A comma seperated list about description of each column of matrix, example: WT_1,WT_2,WT_3,KO_1,KO_2,KO_3 (must be equal to the number of columns containing expression values in matrix)"),
    make_option(c("-w", "--replicates_wt"), default=3, help="Number of replicates for WT treatment (default: %default)"),
    make_option(c("-k", "--replicates_ko"), default=3, help="Number of replicates for KO treatment (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$matrixFile) | is.null(opt$treatment)) {
	cat("\nProgram: matrix2heatmap (R script to plot heatmap corresponding to a heatmap)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

if(identical(opt$matrixFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$matrixFile)
}

treatment <- unlist(strsplit(opt$treatment, ","))
start=2
end=(as.numeric(opt$replicates_wt)+1)
data$WT <- apply(data[,c(start:end)], 1, function(x) mean(x[1:3]))
start=end+1
end=end+as.numeric(opt$replicates_ko)
data$KO <- apply(data[,c(start:end)], 1, function(x) mean(x[1:3]))
data$fc <- log2(data$KO/data$WT)
data <- data[order(data$fc),]
mat <- as.matrix(data[,c(2:end)])
row.names(mat) <- data$V1
colnames(mat) <- treatment
pdf(sprintf("%s.pdf", opt$matrixFile))
heatmap(mat, Rowv=NA, col=cm.colors(256), margins=c(5,10))
dev.off()
