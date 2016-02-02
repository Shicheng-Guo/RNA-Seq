#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="input file containing output from bed2gc_content script")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input))) {
    cat("\nProgram: bed2gc_content.R (plot GC content corresponding to input FASTA sequences\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

data <- read.table(pipe(paste("grep -v ID", opt$input, sep=" ")))
data$per_g <- data$V4/data$V3
data$per_c <- data$V5/data$V3
data$per_a <- data$V6/data$V3
data$per_t <- data$V7/data$V3
data.melt <- melt(data[,c(8:11)])
outfile=sprintf("%s.pdf", opt$input)
g <- ggplot(data.melt, aes(variable, value)) + geom_boxplot(aes(fill=variable)) +
ylim(c(0,1)) + theme_bw()
ggsave(g, dpi=300, filename=outfile, useDingbats=FALSE)

