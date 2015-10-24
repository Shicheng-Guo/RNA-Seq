#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("seqLogo"))
suppressPackageStartupMessages(library("Biostrings"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="input sequences is FASTA format")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input))) {
    cat("\nProgram: fasta2seqlogo.R (plot sequence logo corresponding to input FASTA sequences\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

## credits for the code goes to https://gist.github.com/dianalow/9c9a3b1beed3367300d5
data <- read.table(opt$input)
fasta <- as.vector(data[!grepl(">", data$V1),])
fasta <- fasta[!grepl("N", fasta)]

#plot seq logo
plot_seqlogo <-function(fasta_string){
    require(seqLogo)
    require(Biostrings)
    freq<-consensusMatrix(fasta_string,as.prob=T)[1:4,]
    freq<-data.frame(freq)
    seqLogo(makePWM(freq),ic.scale=FALSE) #ic.scale determines either frequency or bits
    #seqLogo(makePWM(freq),ic.scale=TRUE) #ic.scale determines either frequency or bits
}

outfile=sprintf("%s.pdf", opt$input)
pdf(outfile, width=20)
plot_seqlogo(fasta)
dev.off()

