#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--matrixFile"), help="input count matrix file"),
	make_option(c("-o", "--pdfFile"), help="output pdf file"),
	make_option(c("-s", "--sessionFile"), help="output session file"),
    make_option(c("-n", "--clusterCount"), help="number of clusters")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$matrixFile) | is.null(opt$pdfFile) | is.null(opt$sessionFile) | is.null(opt$clusterCount)) {
	cat("\nProgram: matrix2clustering.R (R script to perform clustering on input count matrix)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(dtw))

if(identical(opt$matrixFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$matrixFile)
}

mat <- as.matrix(data)
set.seed(1)
#kc_script <- kmeans(as.dist(1-cor(t(mat))), as.numeric(opt$clusterCount))
#kc_script <- kmeans(mat, as.numeric(opt$clusterCount))
kc_script <- kmeans(dist(mat, method="DTW"), as.numeric(opt$clusterCount))

mydata <- mat
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)

df <- as.data.frame(mat)
df$clusters <- kc_script$cluster
df <- df[order(df$clusters),]

pdf(opt$pdfFile)
heatmap.2(as.matrix(df[,c(1:ncol(df)-1)]), dendrogram = "none", Rowv = FALSE, Colv = FALSE, trace = "none", col=bluered(10), labRow=NA, scale="row", na.rm=T)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()

save.session(opt$sessionFile)
