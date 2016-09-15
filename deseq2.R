#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--countFile"), help="input read count file (format: <id> <int> <int> <int> (..<int>..) | can be stdin)"),
	make_option(c("-o", "--outDir"), default=".", help="output directory to keep results (default: %default)"),
    make_option(c("-p", "--pvalue"), default="0.05", help="p-value (default: %default)"),
    make_option(c("-n", "--normalized"), help="input file contains normalized read couts", action="store_true"),
    make_option(c("-t", "--treatment"), help="A comma seperated list about description of the data, example: WT,WT,WT,KO,KO,KO (must be equal to the number of samples for which read count is measured in read count file)"),
	make_option(c("-c", "--condition"), help="A comma seperated list about additonal description of the data, example: T1,T2,T2,T1,T2,T2 (must be equal to the number of samples for which read count is measured in read count file)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$countFile) | is.null(opt$treatment)) {
	cat("\nProgram: deseq2.R (R script to identify differential expression using DESeq2)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

## start analysis
## next three lines to debug
#setwd("/Users/sachinpundhir/Lobster/project/rna-seq-analysis/matilda/result_CFUE_2Rep/06_deseq")
#countTable <- read.table("ensembl_transcripts_raw_count", header=TRUE, row.names=1)
#treatment <-unlist(strsplit("WT,WT,WT,KO,KO,KO", ","))
if(identical(opt$countFile, "stdin")==T) {
    countTable <- read.table(file("stdin"), header=F, row.names=1)
} else {
    countTable <- read.table(pipe(paste("grep -v \"^\\_\"", opt$countFile, sep=" ")), header=TRUE, row.names=1)
}

## check number of treatments and conditions, if they are correctly provided
if(ncol(countTable)!=length(unlist(strsplit(opt$treatment, ","))) | (!is.null(opt$condition) & ncol(countTable)!=length(unlist(strsplit(opt$treatment, ","))))) {
  cat("\nProgram: deseq2.R (R script to identify differential expression using DESeq2)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

treatment <- unlist(strsplit(opt$treatment, ","))
colTable <- as.data.frame(as.factor(treatment), row.names=colnames(countTable))
colnames(colTable) <- "treatment"
if(!is.null(opt$condition)) {
  condition <- unlist(strsplit(opt$condition, ","))
  colTable$condition <- as.factor(condition)
}

## check if differential expression is for read count which are a) not normalized, b) normalized
if(is.null(opt$normalized)) {
    ## modify this line to compare more than one set of factors (for example, treatment and condition both)
    if(!is.null(opt$condition)) {
        dds <- DESeqDataSetFromMatrix(countData = countTable, colData = colTable, design = ~ condition + treatment)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = countTable, colData = colTable, design = ~ treatment)
    }
    colData(dds)$treatment <- factor(colData(dds)$treatment, levels=c(levels(colTable$treatment)))
    dds <- DESeq(dds)
    #dds <- DESeq(dds, minReplicatesForReplace=Inf)
} else {
    ## check, if input file is from a cufflinks output
    ## if yes, select only deconvoluted transcripts
    cuff_dir="./";
    if(grepl("/", opt$countFile)) {
        cuff_dir <- gsub("/[^/]+$", "", opt$countFile)
    }

    if(grepl("isoforms", opt$countFile) & file.exists(sprintf("%s/isoforms.count_tracking", cuff_dir))){
        cuff_tracking <- read.table(sprintf("%s/isoforms.count_tracking", cuff_dir), header=TRUE)
        countTable$WT_status <- cuff_tracking$WT_status
        countTable$KO_status <- cuff_tracking$KO_status
        countTable <- countTable[(which(countTable$KO_status=="OK" & countTable$WT_status=="OK")),c(1:6)]
    }

    ## modify this line to compare more than one set of factors (for example, treatment and condition both)
    if(!is.null(opt$condition)) {
        dds <- DESeqDataSetFromMatrix(countData = round(countTable, digits=0), colData = colTable, design = ~ condition + treatment)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = round(countTable, digits=0), colData = colTable, design = ~ treatment)
    }
    colData(dds)$treatment <- factor(colData(dds)$treatment, levels=c(levels(colTable$treatment)))
    normFactors <- matrix(runif(nrow(dds)*ncol(dds),1,1), ncol=ncol(dds),nrow=nrow(dds))
    normalizationFactors(dds) <- normFactors
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
}

## retrieve results
#it makes no difference except that the results() function with no
#arguments will automatically extract results for the last variable in
#the design formula, and for the last level of this variable over the first
#level if this variable is a factor. But otherwise, no it doesn't
#make a difference if you are using the arguments of results() to specify
#which results tables to construct.

## modify this line to compare different set of treatments
res <- results(dds, contrast=c("treatment", "KO", "WT"))
#res <- results(dds, contrast=c("treatment", "t1", "t2"))
#res <- results(dds, contrast=c("treatment", "KO", "WT"), cooksCutoff=FALSE, independentFiltering=FALSE)

## merge results with normalized read counts
size_factors <- sizeFactors(dds)
countTable_norm <- t(apply(countTable, 1, function(x) x/size_factors))
resWithCounts <- merge(res, countTable_norm, by=0, all=TRUE)

resSig <- resWithCounts[which(resWithCounts$padj<as.numeric(opt$pvalue)),]

## define output file name
outFile <- unlist(strsplit(opt$countFile, "/"))[length(unlist(strsplit(opt$countFile, "/")))]

## print results sort by log2foldchange and p-value
## all genes
write.table(as.data.frame(resWithCounts[order(-abs(resWithCounts$log2FoldChange), resWithCounts$padj),]), file=sprintf("%s/%s.all.xls", opt$outDir, outFile), quote=F, sep="\t", row.names=F)

## differentially expressed genes
write.table(as.data.frame(resSig[order(-abs(resSig$log2FoldChange), resSig$padj),]), file=sprintf("%s/%s.de.xls", opt$outDir, outFile), quote=F, sep="\t", row.names=F)

## kind of log output
#print(mcols(res, use.names=TRUE))

## perform quality check
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
hmcol <- colorRampPalette(brewer.pal(8, "RdBu"))(100)
rld <- rlogTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(rownames(colData(dds)), sep=":"))
pdf(sprintf("%s/%s.heatmap.pdf", opt$outDir, outFile))
heatmap.2(mat, trace="none", margin=c(13, 13), col=hmcol)
dev.off()
pdf(sprintf("%s/%s.quality.pdf", opt$outDir, outFile))
par(mfrow=c(2,2))
plotMA(dds, ylim=c(-2,2))
plotDispEsts(dds)
hist(res$pvalue, breaks=20, col="grey")
#plot(attr(res, "filterNumRej"), type="b", xlab="quantile of 'baseMean'", ylab="number of rejections")
dev.off()

#suppressPackageStartupMessages(library(limma)) ## plotMA is available in both DESeq2 and limma
pdf(sprintf("%s/%s.pca.pdf", opt$outDir, outFile))
#par(mfrow=c(2,2))
#plotMDS(log2(assays(dds)[["counts"]]))
#plotMDS(log2(assays(dds)[["mu"]]))
if(!is.null(opt$condition)) {
    data <- plotPCA(rld, intgroup=c("treatment", "condition"), returnData=T)
    percentVar <- round(100 * attr(data, "percentVar"))
    ggplot(data, aes(PC1, PC2, color=treatment, shape=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
} else {
    data <- plotPCA(rld, intgroup="treatment", returnData=T)
    percentVar <- round(100 * attr(data, "percentVar"))
    ggplot(data, aes(PC1, PC2, color=treatment)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
}

dev.off()

## plot volcano plot
pdf(sprintf("%s/%s.volcano.pdf", opt$outDir, outFile), height=15, width=15)
df <- as.data.frame(res)
df$gene <- rownames(df)
df$gene <- gsub(".*_[0-9]+_", "", df$gene)
df <- df[,c(7,2,5,6)]

with(df, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange if log2FC>1, green if both)
#pvalue_threshold=as.numeric(opt$pvalue)
#pvalue_threshold=summary(-log10(res[which(res$padj<0.05),]$pvalue))[5]
pvalue_threshold=1
pvalue_threshold=1/exp(pvalue_threshold*log(10))
#with(subset(df, padj<pvalue_threshold ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#with(subset(df, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(df, padj<pvalue_threshold), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(df, padj<pvalue_threshold), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.8))
dev.off()


## print normalized read count of genes
## this was WRONG (https://support.bioconductor.org/p/66067/)
#write.table(as.data.frame(assay(rld)), file=sprintf("%s/%s.normExpr.xls", opt$outDir, outFile), sep="\t", quote=F, col.names=T, row.names=T)
## This is RIGHT (use following instead) (sorted by log2fold change)
write.table(countTable_norm[order(-as.data.frame(res)[,2]),], file=sprintf("%s/%s.normExpr.xls", opt$outDir, outFile), sep="\t", quote=F, col.names=F, row.names=T)

## print size factors for samples
print(size_factors)

## save the session for further use
save.session(sprintf("%s/%s.session", opt$outDir, outFile))


## code for unknown batch effects
#If you have known batches, just include the batch variable in the design for DESeq2. 

#We don't recommend testing on transformed counts.

#If you have unknown batches, you can use svaseq or other packages. We are writing up a workflow which will be released in a few weeks and includes svaseq and RUVSeq. 

#But briefly, add the SVA surrogate variables (columns of 'sv') to the colData, and then add these to the design. E.g., for two surrogate variables:

#Code:
#dds$SV1 <- svseq$sv[,1]
#dds$SV2 <- svseq$sv[,2]
#design(dds) <- ~ SV1 + SV2 + condition
#dds <- DESeq(dds)
## normalized counts after batch correction
#fitted.common.scale = t( t( assays(dds)[["mu"]] ) / sizeFactors(dds) )
