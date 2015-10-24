#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--id"), help="any user defined unique id"),
	make_option(c("-r", "--bgFile"), help="block group statistics file"),
	make_option(c("-b", "--bFile"), help="block statistics file"),
	make_option(c("-o", "--outDir"), help="output directory", default=".")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$id) | is.null(opt$bgFile) | is.null(opt$bFile)) {
	cat("\nProgram: plotBlockStat.r (plot block group and block statistics in pdf format)\n")
	cat("Author: University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: sachin@rth.dk\n")
	print_help(parser)
	q()
}

bg <- read.table(opt$bgFile, sep="\t")
b <- read.table(opt$bFile, sep="\t")

## compute block group size
if(is.null(bg$V10)) {
	colFirst <- gsub("^.+chr", "", bg[,1]);
	size <- length(unlist(strsplit(as.character(colFirst[1]), "_", fixed=TRUE)));
	tmp <- unlist(strsplit(as.character(colFirst), "_", fixed=TRUE));
	matrix <- matrix(tmp, nrow=length(tmp)/size, ncol=size, byrow=TRUE);
	bg$V10 <- (as.numeric(matrix[,(size-1)]) - as.numeric(matrix[,(size-2)]))+1;
}
## compute percentage unique reads
bg$V11 <- (bg$V9/bg$V3)*100
## read annotation for block group
bg.mirna <- bg[which(bg$V5=="miRNA"),]
bg.snorna <- bg[grep("snoRNA*", bg$V5),]
bg.trna <- bg[which(bg$V5=="tRNA"),]
bg.snrna <- bg[which(bg$V5=="snRNA"),]
bg.scrna <- bg[which(bg$V5=="scRNA"),]
bg.yrna <- bg[which(bg$V5=="Y_RNA"),]
bg.na <- bg[which(bg$V5=="n/a"),]
bg.others <- bg[setdiff(seq(length(bg$V5)), grep("snoRNA*", bg$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("miRNA", bg.others$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("tRNA", bg.others$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("snRNA", bg.others$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("scRNA", bg.others$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("Y_RNA", bg.others$V5)),]
bg.others <- bg.others[setdiff(seq(length(bg.others$V5)), grep("n/a", bg.others$V5)),]
## read loci for block group
bg.exon <- bg[grep("exon.*", bg$V6),]
bg.intron <- bg[grep("intron.*", bg$V6),]
bg.utr3 <- bg[grep("UTR3.*", bg$V6),]
bg.utr5 <- bg[grep("UTR5.*", bg$V6),]
bg.intergene <- bg[setdiff(seq(length(bg$V6)), grep("exon.*", bg$V6)),]
bg.intergene <- bg.intergene[setdiff(seq(length(bg.intergene$V6)), grep("intron.*", bg.intergene$V6)),]
bg.intergene <- bg.intergene[setdiff(seq(length(bg.intergene$V6)), grep("UTR5.*", bg.intergene$V6)),]
bg.intergene <- bg.intergene[setdiff(seq(length(bg.intergene$V6)), grep("UTR3.*", bg.intergene$V6)),]
## read annotation for block
b.mirna <- b[which(b$V7=="miRNA"),]
b.snorna <- b[grep("snoRNA*", b$V7),]
b.trna <- b[which(b$V7=="tRNA"),]
b.snrna <- b[which(b$V7=="snRNA"),]
b.scrna <- b[which(b$V7=="scRNA"),]
b.yrna <- b[which(b$V7=="Y_RNA"),]
b.na <- b[which(b$V7=="n/a"),]
b.others <- b[setdiff(seq(length(b$V7)), grep("snoRNA*", b$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("miRNA", b.others$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("tRNA", b.others$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("snRNA", b.others$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("scRNA", b.others$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("Y_RNA", b.others$V7)),]
b.others <- b.others[setdiff(seq(length(b.others$V7)), grep("n/a", b.others$V7)),]
## read loci for block
b.exon <- b[grep("exon.*", b$V8),]
b.intron <- b[grep("intron.*", b$V8),]
b.utr3 <- b[grep("UTR3.*", b$V8),]
b.utr5 <- b[grep("UTR5.*", b$V8),]
b.intergene <- b[setdiff(seq(length(b$V8)), grep("exon.*", b$V8)),]
b.intergene <- b.intergene[setdiff(seq(length(b.intergene$V8)), grep("intron.*", b.intergene$V8)),]
b.intergene <- b.intergene[setdiff(seq(length(b.intergene$V8)), grep("UTR5.*", b.intergene$V8)),]
b.intergene <- b.intergene[setdiff(seq(length(b.intergene$V8)), grep("UTR3.*", b.intergene$V8)),]


infile <- paste(opt$outDir, "/", opt$id, "Bar.eps", sep="");
#postscript(file=infile, paper="special", width=7, height=6.5, horizontal=F)
postscript(file=infile, horizontal=F)
par(mfrow=c(2,1))
# barplot: block count
breaks <- seq(1,10, by=1)
bg.mirna.density <- c(length(bg.mirna[which(bg.mirna$V2==1),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==2),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==3),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==4),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==5),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==6),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==7),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==8),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==9),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2==10),]$V2)/nrow(bg.mirna), length(bg.mirna[which(bg.mirna$V2>10),]$V2)/nrow(bg.mirna))
bg.snorna.density <- c(length(bg.snorna[which(bg.snorna$V2==1),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==2),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==3),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==4),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==5),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==6),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==7),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==8),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==9),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2==10),]$V2)/nrow(bg.snorna), length(bg.snorna[which(bg.snorna$V2>10),]$V2)/nrow(bg.snorna))
bg.trna.density <- c(length(bg.trna[which(bg.trna$V2==1),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==2),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==3),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==4),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==5),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==6),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==7),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==8),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==9),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2==10),]$V2)/nrow(bg.trna), length(bg.trna[which(bg.trna$V2>10),]$V2)/nrow(bg.trna))
bg.snrna.density <- c(length(bg.snrna[which(bg.snrna$V2==1),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==2),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==3),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==4),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==5),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==6),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==7),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==8),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==9),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2==10),]$V2)/nrow(bg.snrna), length(bg.snrna[which(bg.snrna$V2>10),]$V2)/nrow(bg.snrna))
bg.scrna.density <- c(length(bg.scrna[which(bg.scrna$V2==1),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==2),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==3),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==4),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==5),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==6),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==7),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==8),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==9),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2==10),]$V2)/nrow(bg.scrna), length(bg.scrna[which(bg.scrna$V2>10),]$V2)/nrow(bg.scrna))
bg.yrna.density <- c(length(bg.yrna[which(bg.yrna$V2==1),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==2),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==3),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==4),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==5),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==6),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==7),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==8),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==9),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2==10),]$V2)/nrow(bg.yrna), length(bg.yrna[which(bg.yrna$V2>10),]$V2)/nrow(bg.yrna))
bg.others.density <- c(length(bg.others[which(bg.others$V2==1),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==2),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==3),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==4),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==5),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==6),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==7),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==8),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==9),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2==10),]$V2)/nrow(bg.others), length(bg.others[which(bg.others$V2>10),]$V2)/nrow(bg.others))
bg.na.density <- c(length(bg.na[which(bg.na$V2==1),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==2),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==3),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==4),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==5),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==6),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==7),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==8),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==9),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2==10),]$V2)/nrow(bg.na), length(bg.na[which(bg.na$V2>10),]$V2)/nrow(bg.na))
barplot(rbind(bg.mirna.density, bg.snorna.density, bg.trna.density, bg.snrna.density, bg.scrna.density, bg.others.density, bg.na.density, bg.yrna.density), beside=1, names.arg=as.character(c(breaks[1:10], ">10")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="Block count", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)

## barplot: unique reads percentage
breaks <- seq(0,100,by=10)
bg.mirna.freq <- hist(bg.mirna$V11, breaks, right=FALSE, plot=FALSE)
bg.snorna.freq <- hist(bg.snorna$V11, breaks, right=FALSE, plot=FALSE)
bg.trna.freq <- hist(bg.trna$V11, breaks, right=FALSE, plot=FALSE)
bg.snrna.freq <- hist(bg.snrna$V11, breaks, right=FALSE, plot=FALSE)
bg.scrna.freq <- hist(bg.scrna$V11, breaks, right=FALSE, plot=FALSE)
bg.yrna.freq <- hist(bg.yrna$V11, breaks, right=FALSE, plot=FALSE)
bg.others.freq <- hist(bg.others$V11, breaks, right=FALSE, plot=FALSE)
bg.na.freq <- hist(bg.na$V11, breaks, right=FALSE, plot=FALSE)
barplot(rbind(bg.mirna.freq$density*10, bg.snorna.freq$density*10, bg.trna.freq$density*10, bg.snrna.freq$density*10, bg.scrna.freq$density*10, bg.na.freq$density*10, bg.others.freq$density*10, bg.yrna.freq$density*10), beside=1, names.arg=as.character(c(breaks[2:length(breaks)])), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="unique reads (%)", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)

par(mfrow=c(2,1))
## barplot: block group entropy
breaks <- seq(0,round(max(bg$V4))+1, by=1)
bg.mirna.freq <- hist(bg.mirna$V4, breaks, right=FALSE, plot=FALSE)
bg.snorna.freq <- hist(bg.snorna$V4, breaks, right=FALSE, plot=FALSE)
bg.trna.freq <- hist(bg.trna$V4, breaks, right=FALSE, plot=FALSE)
bg.snrna.freq <- hist(bg.snrna$V4, breaks, right=FALSE, plot=FALSE)
bg.scrna.freq <- hist(bg.scrna$V4, breaks, right=FALSE, plot=FALSE)
bg.yrna.freq <- hist(bg.yrna$V4, breaks, right=FALSE, plot=FALSE)
bg.others.freq <- hist(bg.others$V4, breaks, right=FALSE, plot=FALSE)
bg.na.freq <- hist(bg.na$V4, breaks, right=FALSE, plot=FALSE)
barplot(rbind(bg.mirna.freq$density, bg.snorna.freq$density, bg.trna.freq$density, bg.snrna.freq$density, bg.scrna.freq$density, bg.na.freq$density, bg.others.freq$density, bg.yrna.freq$density), beside=1, names.arg=as.character(c(breaks[2:length(breaks)])), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="Entropy", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)

## barplot: block group relEntropy
breaks <- seq(0,round(max(bg$V8))+1, by=1)
bg.mirna.freq <- hist(bg.mirna$V8, breaks, right=FALSE, plot=FALSE)
bg.snorna.freq <- hist(bg.snorna$V8, breaks, right=FALSE, plot=FALSE)
bg.trna.freq <- hist(bg.trna$V8, breaks, right=FALSE, plot=FALSE)
bg.snrna.freq <- hist(bg.snrna$V8, breaks, right=FALSE, plot=FALSE)
bg.scrna.freq <- hist(bg.scrna$V8, breaks, right=FALSE, plot=FALSE)
bg.yrna.freq <- hist(bg.yrna$V8, breaks, right=FALSE, plot=FALSE)
bg.others.freq <- hist(bg.others$V8, breaks, right=FALSE, plot=FALSE)
bg.na.freq <- hist(bg.na$V8, breaks, right=FALSE, plot=FALSE)
barplot(rbind(bg.mirna.freq$density, bg.snorna.freq$density, bg.trna.freq$density, bg.snrna.freq$density, bg.scrna.freq$density, bg.na.freq$density, bg.others.freq$density, bg.yrna.freq$density), beside=1, names.arg=as.character(c(breaks[2:length(breaks)])), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="Relative entropy", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)

par(mfrow=c(2,1))
## barplot: block group loci
slices <- c(nrow(bg.mirna[grep("exon.*", bg.mirna$V6),]), nrow(bg.mirna[grep("intron.*", bg.mirna$V6),]), nrow(bg.mirna[grep("UTR5.*", bg.mirna$V6),]), nrow(bg.mirna[grep("UTR3.*", bg.mirna$V6),]), (nrow(bg.mirna)-(nrow(bg.mirna[grep("exon.*", bg.mirna$V6),])+nrow(bg.mirna[grep("intron.*", bg.mirna$V6),])+nrow(bg.mirna[grep("UTR5.*", bg.mirna$V6),])+nrow(bg.mirna[grep("UTR3.*", bg.mirna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.mirna <- slices/sum(slices)
lbls <- paste(lbls, pct.mirna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="a) miRNA")
slices <- c(nrow(bg.snorna[grep("exon.*", bg.snorna$V6),]), nrow(bg.snorna[grep("intron.*", bg.snorna$V6),]), nrow(bg.snorna[grep("UTR5.*", bg.snorna$V6),]), nrow(bg.snorna[grep("UTR3.*", bg.snorna$V6),]), (nrow(bg.snorna)-(nrow(bg.snorna[grep("exon.*", bg.snorna$V6),])+nrow(bg.snorna[grep("intron.*", bg.snorna$V6),])+nrow(bg.snorna[grep("UTR5.*", bg.snorna$V6),])+nrow(bg.snorna[grep("UTR3.*", bg.snorna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.snorna <- slices/sum(slices)
lbls <- paste(lbls, pct.snorna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="b) snoRNA")
slices <- c(nrow(bg.trna[grep("exon.*", bg.trna$V6),]), nrow(bg.trna[grep("intron.*", bg.trna$V6),]), nrow(bg.trna[grep("UTR5.*", bg.trna$V6),]), nrow(bg.trna[grep("UTR3.*", bg.trna$V6),]), (nrow(bg.trna)-(nrow(bg.trna[grep("exon.*", bg.trna$V6),])+nrow(bg.trna[grep("intron.*", bg.trna$V6),])+nrow(bg.trna[grep("UTR5.*", bg.trna$V6),])+nrow(bg.trna[grep("UTR3.*", bg.trna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.trna <- slices/sum(slices)
lbls <- paste(lbls, pct.trna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="c) tRNA")
slices <- c(nrow(bg.snrna[grep("exon.*", bg.snrna$V6),]), nrow(bg.snrna[grep("intron.*", bg.snrna$V6),]), nrow(bg.snrna[grep("UTR5.*", bg.snrna$V6),]), nrow(bg.snrna[grep("UTR3.*", bg.snrna$V6),]), (nrow(bg.snrna)-(nrow(bg.snrna[grep("exon.*", bg.snrna$V6),])+nrow(bg.snrna[grep("intron.*", bg.snrna$V6),])+nrow(bg.snrna[grep("UTR5.*", bg.snrna$V6),])+nrow(bg.snrna[grep("UTR3.*", bg.snrna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.snrna <- slices/sum(slices)
lbls <- paste(lbls, pct.snrna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="d) snRNA")
slices <- c(nrow(bg.scrna[grep("exon.*", bg.scrna$V6),]), nrow(bg.scrna[grep("intron.*", bg.scrna$V6),]), nrow(bg.scrna[grep("UTR5.*", bg.scrna$V6),]), nrow(bg.scrna[grep("UTR3.*", bg.scrna$V6),]), (nrow(bg.scrna)-(nrow(bg.scrna[grep("exon.*", bg.scrna$V6),])+nrow(bg.scrna[grep("intron.*", bg.scrna$V6),])+nrow(bg.scrna[grep("UTR5.*", bg.scrna$V6),])+nrow(bg.scrna[grep("UTR3.*", bg.scrna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.scrna <- slices/sum(slices)
lbls <- paste(lbls, pct.scrna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="e) scRNA")
slices <- c(nrow(bg.yrna[grep("exon.*", bg.yrna$V6),]), nrow(bg.yrna[grep("intron.*", bg.yrna$V6),]), nrow(bg.yrna[grep("UTR5.*", bg.yrna$V6),]), nrow(bg.yrna[grep("UTR3.*", bg.yrna$V6),]), (nrow(bg.yrna)-(nrow(bg.yrna[grep("exon.*", bg.yrna$V6),])+nrow(bg.yrna[grep("intron.*", bg.yrna$V6),])+nrow(bg.yrna[grep("UTR5.*", bg.yrna$V6),])+nrow(bg.yrna[grep("UTR3.*", bg.yrna$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.yrna <- slices/sum(slices)
lbls <- paste(lbls, pct.yrna)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="f) yRNA")
slices <- c(nrow(bg.others[grep("exon.*", bg.others$V6),]), nrow(bg.others[grep("intron.*", bg.others$V6),]), nrow(bg.others[grep("UTR5.*", bg.others$V6),]), nrow(bg.others[grep("UTR3.*", bg.others$V6),]), (nrow(bg.others)-(nrow(bg.others[grep("exon.*", bg.others$V6),])+nrow(bg.others[grep("intron.*", bg.others$V6),])+nrow(bg.others[grep("UTR5.*", bg.others$V6),])+nrow(bg.others[grep("UTR3.*", bg.others$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.others <- slices/sum(slices)
lbls <- paste(lbls, pct.others)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="g) other ncRNAs")
slices <- c(nrow(bg.na[grep("exon.*", bg.na$V6),]), nrow(bg.na[grep("intron.*", bg.na$V6),]), nrow(bg.na[grep("UTR5.*", bg.na$V6),]), nrow(bg.na[grep("UTR3.*", bg.na$V6),]), (nrow(bg.na)-(nrow(bg.na[grep("exon.*", bg.na$V6),])+nrow(bg.na[grep("intron.*", bg.na$V6),])+nrow(bg.na[grep("UTR5.*", bg.na$V6),])+nrow(bg.na[grep("UTR3.*", bg.na$V6),]))))
lbls <- c("Exon", "Intron", "UTR5", "UTR3", "Intergene")
pct.na <- slices/sum(slices)
lbls <- paste(lbls, pct.na)
lbls <- paste(lbls, "%", sep="")
#pie(slices, lbls, main="h) unannotated block groups")
pct.exon <- c(pct.mirna[1], pct.snorna[1], pct.trna[1], pct.snrna[1], pct.scrna[1], pct.yrna[1], pct.others[1], pct.na[1])
pct.intron <- c(pct.mirna[2], pct.snorna[2], pct.trna[2], pct.snrna[2], pct.scrna[2], pct.yrna[2], pct.others[2], pct.na[2])
pct.utr5 <- c(pct.mirna[3], pct.snorna[3], pct.trna[3], pct.snrna[3], pct.scrna[3], pct.yrna[3], pct.others[3], pct.na[3])
pct.utr3 <- c(pct.mirna[4], pct.snorna[4], pct.trna[4], pct.snrna[4], pct.scrna[4], pct.yrna[4], pct.others[4], pct.na[4])
pct.intergene <- c(pct.mirna[5], pct.snorna[5], pct.trna[5], pct.snrna[5], pct.scrna[5], pct.yrna[5], pct.others[5], pct.na[5])
barplot(rbind(pct.exon, pct.intron, pct.utr5, pct.utr3, pct.intergene), beside=1, names.arg=as.character(c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3"), xlab="Annotation", ylab="Density", ylim=c(0,1), cex.axis=1.0, cex.lab=1.0, cex.names=1.0)
legend("topright", legend=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3"), cex=1.2)

## plot to show relation between entropy, block group size and block count
bgCount1 <- bg[which(bg$V2<=2 & bg$V10<=500),]
bgCount2 <- bg[which(bg$V2>2 & bg$V2 <=4 & bg$V10<=500),]
bgCount3 <- bg[which(bg$V2>4 & bg$V10<=500),]
plot(bg[which(bg$V10<=500),]$V10, bg[which(bg$V10<=500),]$V4, type="n", xlab="block group size", ylab="entropy")
points(bgCount1$V10, bgCount1$V4, pch=1, col="goldenrod4")
points(bgCount2$V10, bgCount2$V4, pch=1, col="lightgray")
points(bgCount3$V10, bgCount3$V4, pch=1, col="coral")
legend("topright", title="block count", legend=c("<=2", "3-4 ", ">4"), fill=c("goldenrod4", "lightgray", "coral"), cex=1.2)

par(mfrow=c(2,1))
## barplot: block size
breaks <- seq(0,300,by=50)
b.mirna.freq <- hist(b.mirna[which(b.mirna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.mirna.freq$density <- c(b.mirna.freq$counts/length(b.mirna$V2), length(b.mirna[which(b.mirna$V2>300),]$V2)/length(b.mirna$V2))
b.snorna.freq <- hist(b.snorna[which(b.snorna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.snorna.freq$density <- c(b.snorna.freq$counts/length(b.snorna$V2), length(b.snorna[which(b.snorna$V2>300),]$V2)/length(b.snorna$V2))
b.trna.freq <- hist(b.trna[which(b.trna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.trna.freq$density <- c(b.trna.freq$counts/length(b.trna$V2), length(b.trna[which(b.trna$V2>300),]$V2)/length(b.trna$V2))
b.snrna.freq <- hist(b.snrna[which(b.snrna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.snrna.freq$density <- c(b.snrna.freq$counts/length(b.snrna$V2), length(b.snrna[which(b.snrna$V2>300),]$V2)/length(b.snrna$V2))
b.scrna.freq <- hist(b.scrna[which(b.scrna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.scrna.freq$density <- c(b.scrna.freq$counts/length(b.scrna$V2), length(b.scrna[which(b.scrna$V2>300),]$V2)/length(b.scrna$V2))
b.yrna.freq <- hist(b.yrna[which(b.yrna$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.yrna.freq$density <- c(b.yrna.freq$counts/length(b.yrna$V2), length(b.yrna[which(b.yrna$V2>300),]$V2)/length(b.yrna$V2))
b.others.freq <- hist(b.others[which(b.others$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.others.freq$density <- c(b.others.freq$counts/length(b.others$V2), length(b.others[which(b.others$V2>300),]$V2)/length(b.others$V2))
b.na.freq <- hist(b.na[which(b.na$V2<=300),]$V2, breaks, right=FALSE, plot=FALSE)
b.na.freq$density <- c(b.na.freq$counts/length(b.na$V2), length(b.na[which(b.na$V2>300),]$V2)/length(b.na$V2))
barplot(rbind(b.mirna.freq$density, b.snorna.freq$density, b.trna.freq$density, b.snrna.freq$density, b.scrna.freq$density, b.na.freq$density, b.others.freq$density, b.yrna.freq$density), beside=1, names.arg=as.character(c(breaks[2:length(breaks)], ">300")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="Block size", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)

## barplot: block group size
breaks <- seq(0,500,by=50)
#j=1; for(i in (bg[,1])) { tmp <- unlist(strsplit(as.character(i), "_", fixed=TRUE)); bg[j,9] <- (as.numeric(tmp[3])-as.numeric(tmp[2]))+1; j=j+1; }
bg.mirna.freq <- hist(bg.mirna[which(bg.mirna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.mirna.freq$density <- c(bg.mirna.freq$counts/length(bg.mirna$V10), length(bg.mirna[which(bg.mirna$V10>500),]$V10)/length(bg.mirna$V10))
bg.mirna.freq$counts <- c(bg.mirna.freq$counts, length(bg.mirna[which(bg.mirna$V10>500),]$V10))
bg.snorna.freq <- hist(bg.snorna[which(bg.snorna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.snorna.freq$density <- c(bg.snorna.freq$counts/length(bg.snorna$V10), length(bg.snorna[which(bg.snorna$V10>500),]$V10)/length(bg.snorna$V10))
bg.snorna.freq$counts <- c(bg.snorna.freq$counts, length(bg.snorna[which(bg.snorna$V10>500),]$V10))
bg.trna.freq <- hist(bg.trna[which(bg.trna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.trna.freq$density <- c(bg.trna.freq$counts/length(bg.trna$V10), length(bg.trna[which(bg.trna$V10>500),]$V10)/length(bg.trna$V10))
bg.trna.freq$counts <- c(bg.trna.freq$counts, length(bg.trna[which(bg.trna$V10>500),]$V10))
bg.snrna.freq <- hist(bg.snrna[which(bg.snrna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.snrna.freq$density <- c(bg.snrna.freq$counts/length(bg.snrna$V10), length(bg.snrna[which(bg.snrna$V10>500),]$V10)/length(bg.snrna$V10))
bg.snrna.freq$counts <- c(bg.snrna.freq$counts, length(bg.snrna[which(bg.snrna$V10>500),]$V10))
bg.scrna.freq <- hist(bg.scrna[which(bg.scrna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.scrna.freq$density <- c(bg.scrna.freq$counts/length(bg.scrna$V10), length(bg.scrna[which(bg.scrna$V10>500),]$V10)/length(bg.scrna$V10))
bg.scrna.freq$counts <- c(bg.scrna.freq$counts, length(bg.scrna[which(bg.scrna$V10>500),]$V10))
bg.yrna.freq <- hist(bg.yrna[which(bg.yrna$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.yrna.freq$density <- c(bg.yrna.freq$counts/length(bg.yrna$V10), length(bg.yrna[which(bg.yrna$V10>500),]$V10)/length(bg.yrna$V10))
bg.yrna.freq$counts <- c(bg.yrna.freq$counts, length(bg.yrna[which(bg.yrna$V10>500),]$V10))
bg.others.freq <- hist(bg.others[which(bg.others$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.others.freq$density <- c(bg.others.freq$counts/length(bg.others$V10), length(bg.others[which(bg.others$V10>500),]$V10)/length(bg.others$V10))
bg.others.freq$counts <- c(bg.others.freq$counts, length(bg.others[which(bg.others$V10>500),]$V10))
bg.na.freq <- hist(bg.na[which(bg.na$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.na.freq$density <- c(bg.na.freq$counts/length(bg.na$V10), length(bg.na[which(bg.na$V10>500),]$V10)/length(bg.na$V10))
bg.na.freq$counts <- c(bg.na.freq$counts, length(bg.na[which(bg.na$V10>500),]$V10))
bg.freq <- hist(bg[which(bg$V10<=500),]$V10, breaks, right=FALSE, plot=FALSE)
bg.freq$density <- c(bg.freq$counts/length(bg$V10), length(bg[which(bg$V10>500),]$V10)/length(bg$V10))
bg.freq$counts <- c(bg.freq$counts, length(bg[which(bg$V10>500),]$V10))
barplot(rbind(bg.mirna.freq$density, bg.snorna.freq$density, bg.trna.freq$density, bg.snrna.freq$density, bg.scrna.freq$density, bg.na.freq$density, bg.others.freq$density, bg.yrna.freq$density, bg.freq$density), beside=1, names.arg=as.character(c(breaks[2:length(breaks)], ">500")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray", "black"), xlab="Block group size", ylab="Density", ylim=c(0,1), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
#barplot(rbind(log(bg.mirna.freq$counts), log(bg.snorna.freq$counts), log(bg.trna.freq$counts), log(bg.others.freq$counts), log(bg.na.freq$counts), log(bg.freq$counts)), beside=1, names.arg=as.character(c(breaks[2:length(breaks)], ">500")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "coral"), xlab="Block size", ylab="Density", ylim=c(0,10), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA", "All"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray", "black"), cex=1.2)

par(mfrow=c(2,1))
## barplot block distance
breaks <- seq(0,200,by=20)
b.mirna.freq <- hist(b.mirna[which(b.mirna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.mirna.freq$density <- c(b.mirna.freq$counts/length(which(is.na(b.mirna$V3)==FALSE)), length(b.mirna[which(b.mirna$V3>200),]$V3)/length(which(is.na(b.mirna$V3)==FALSE)));
b.snorna.freq <- hist(b.snorna[which(b.snorna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.snorna.freq$density <- c(b.snorna.freq$counts/length(which(is.na(b.snorna$V3)==FALSE)), length(b.snorna[which(b.snorna$V3>200),]$V3)/length(which(is.na(b.snorna$V3)==FALSE)));
b.trna.freq <- hist(b.trna[which(b.trna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.trna.freq$density <- c(b.trna.freq$counts/length(which(is.na(b.trna$V3)==FALSE)), length(b.trna[which(b.trna$V3>200),]$V3)/length(which(is.na(b.trna$V3)==FALSE)));
b.snrna.freq <- hist(b.snrna[which(b.snrna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.snrna.freq$density <- c(b.snrna.freq$counts/length(which(is.na(b.snrna$V3)==FALSE)), length(b.snrna[which(b.snrna$V3>200),]$V3)/length(which(is.na(b.snrna$V3)==FALSE)));
b.scrna.freq <- hist(b.scrna[which(b.scrna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.scrna.freq$density <- c(b.scrna.freq$counts/length(which(is.na(b.scrna$V3)==FALSE)), length(b.scrna[which(b.scrna$V3>200),]$V3)/length(which(is.na(b.scrna$V3)==FALSE)));
b.yrna.freq <- hist(b.yrna[which(b.yrna$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.yrna.freq$density <- c(b.yrna.freq$counts/length(which(is.na(b.yrna$V3)==FALSE)), length(b.yrna[which(b.yrna$V3>200),]$V3)/length(which(is.na(b.yrna$V3)==FALSE)));
b.others.freq <- hist(b.others[which(b.others$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.others.freq$density <- c(b.others.freq$counts/length(which(is.na(b.others$V3)==FALSE)), length(b.others[which(b.others$V3>200),]$V3)/length(which(is.na(b.others$V3)==FALSE)));
b.na.freq <- hist(b.na[which(b.na$V3<=200),]$V3, breaks, right=FALSE, plot=FALSE)
b.na.freq$density <- c(b.na.freq$counts/length(which(is.na(b.na$V3)==FALSE)), length(b.na[which(b.na$V3>200),]$V3)/length(which(is.na(b.na$V3)==FALSE)));
barplot(rbind(b.mirna.freq$density, b.snorna.freq$density, b.trna.freq$density, b.snrna.freq$density, b.scrna.freq$density, b.na.freq$density, b.others.freq$density, b.yrna.freq$density), beside=1, names.arg=as.character(c(breaks[2:length(breaks)],">200")), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), xlab="Block distance", ylab="Density", ylim=c(0,0.5), cex.axis=1.2, cex.lab=1.2, cex.names=1.2)
legend("topright", legend=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Other", "UA", "Y_RNA"), fill=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), cex=1.2)
dev.off()


infile <- paste(opt$outDir, "/", opt$id, "Stat.eps", sep="");
#pdf(file=infile, width=10, height=20, onefile=TRUE, family='Helvetica', paper='letter', pointsize=8)
postscript(file=infile, horizontal=F)
## plot annotation for block group
par(mfrow=c(2,3))
boxplot(list(bg.mirna$V2, bg.snorna$V2, bg.trna$V2, bg.snrna$V2, bg.scrna$V2, bg.yrna$V2, bg.others$V2, bg.na$V2), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Block Count", font=2))
boxplot(list(bg.mirna$V3, bg.snorna$V3, bg.trna$V3, bg.snrna$V3, bg.scrna$V3, bg.yrna$V3, bg.others$V3, bg.na$V3), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE, log="y")
title(main=list("Total Reads", font=2))
boxplot(list(bg.mirna$V4, bg.snorna$V4, bg.trna$V4, bg.snrna$V4, bg.scrna$V4, bg.yrna$V4, bg.others$V4, bg.na$V4), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Entropy", font=2))
## plot loci for block group
boxplot(list(bg.exon$V2, bg.intron$V2, bg.utr5$V2, bg.utr3$V2, bg.intergene$V2), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Block Count", font=2))
boxplot(list(bg.exon$V3, bg.intron$V3, bg.utr5$V3, bg.utr3$V3, bg.intergene$V3), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Total Reads", font=2))
boxplot(list(bg.exon$V4, bg.intron$V4, bg.utr5$V4, bg.utr3$V4, bg.intergene$V4), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Block Group Entropy", font=2))
## plot annotation for block
par(mfrow=c(2,3))
boxplot(list(b.mirna$V2, b.snorna$V2, b.trna$V2, b.snrna$V2, b.scrna$V2, b.yrna$V2, b.others$V2, b.na$V2), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Block Size", font=2))
boxplot(list(b.mirna$V3, b.snorna$V3, b.trna$V3, b.snrna$V3, b.scrna$V3, b.yrna$V3, b.others$V3, b.na$V3), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Block Distance", font=2))
boxplot(list(b.mirna$V4, b.snorna$V4, b.trna$V4, b.snrna$V4, b.scrna$V4, b.yrna$V4, b.others$V4, b.na$V4), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Block Entropy", font=2))
boxplot(list(b.mirna$V5, b.snorna$V5, b.trna$V5, b.snrna$V5, b.scrna$V5, b.yrna$V5, b.others$V5, b.na$V5), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Total Reads", font=2))
boxplot(list(b.mirna$V6, b.snorna$V6, b.trna$V6, b.snrna$V6, b.scrna$V6, b.yrna$V6, b.others$V6, b.na$V6), names=c("miRNA", "snoRNA", "tRNA", "snRNA", "scRNA", "Y_RNA", "Other", "UA"), col=c("khaki3", "olivedrab4", "lightblue3", "slategray4", "snow3", "mistyrose3", "lightsalmon3", "lightgray"), outline=FALSE)
title(main=list("Overlap", font=2))
## plot loci for block
par(mfrow=c(2,3))
boxplot(list(b.exon$V2, b.intron$V2, b.utr5$V2, b.utr3$V2, b.intergene$V2), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Block Size", font=2))
boxplot(list(b.exon$V3, b.intron$V3, b.utr5$V3, b.utr3$V3, b.intergene$V3), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Block Distance", font=2))
boxplot(list(b.exon$V4, b.intron$V4, b.utr5$V4, b.utr3$V4, b.intergene$V4), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Block Entropy", font=2))
boxplot(list(b.exon$V5, b.intron$V5, b.utr5$V5, b.utr3$V5, b.intergene$V5), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Total Reads", font=2))
boxplot(list(b.exon$V6, b.intron$V6, b.utr5$V6, b.utr3$V6, b.intergene$V6), names=c("Exon", "Intron", "UTR5", "UTR3", "Intergene"), col=c("lightsteelblue3", "lightgray", "mistyrose3", "navajowhite3", "lightsalmon3"), outline=FALSE)
title(main=list("Overlap", font=2))
dev.off()
q()
