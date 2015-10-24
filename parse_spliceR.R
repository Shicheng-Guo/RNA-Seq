#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("spliceR"))
suppressPackageStartupMessages(library("session"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("WriteXLS"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="input spliceR directory having Rdata"),
    make_option(c("-q", "--qvalue"), help="minimum q-value to output genes or transcripts as up- and down-regulated (default: %default)", default=0.05),
    make_option(c("-f", "--foldchange"), help="if provided, log fold change is used in place of q-value to determine up- and down-regulated genes"),
    make_option(c("-x", "--minexpr"), help="minimum expression (FPKM) to consider a gene or transcript as up- or down-regulated (default: %default). Only used if -f is defined. This value is also used to compute transcripts specific to KO or WT samples, irrespective of whether -f is defined or not", default=5),
    make_option(c("-y", "--minisoexpr"), help="minimum expression (FPKM) of isoforms to analyze them for isoform switch (default: %default)", default=30),
    make_option(c("-z", "--minisoexprdiff"), help="minimum difference in the expression (FPKM) between two isoforms to define them as isoform switch (default: %default)", default=10),
    make_option(c("-t", "--generateGTF"), help="generate GTF files also", action="store_true"),
    make_option(c("-m", "--organism"), help="genome (default: %default)", default="mm9"),
    #make_option(c("-d", "--deseq2ResFile"), help="input file containing complete results from DESeq2. If provided, the results will undergo independent filtering during multiple correction"),
    make_option(c("-d", "--indFilter"), help="input directory containing results from cuffdiff. If provided, it will be used to recompute q-values based on Independent filtering. This is similar to as performed in DESeq2"),
    make_option(c("-e", "--deseq2SigGeneFile"), help="input file containing differentially expressed genes from DESeq2. If provided, a venn diagram showing overlap between the two prediction methods (cuffdiff and DESeq2) will be plotted."),
    make_option(c("-s", "--deseq2SigIsoFile"), help="input file containing differentially expressed transcripts from DESeq2. If provided, a venn diagram showing overlap between the two prediction methods (cuffdiff and DESeq2) will be plotted."),
    make_option(c("-l", "--lowConfAlso"), help="include low confidence predictions also (default is to include only high confidence)", action="store_true"),
    make_option(c("-r", "--xlsFormat"), help="output file in XLS format (default: table)", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input))) {
    cat("\nProgram: parse_spliceR (parse spliceR output files to determine up-, down-regulated and PTC genes or isoforms)\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

setwd(opt$input)

# load data
if(is.null(opt$lowConfAlso)) {
    load("03_mySpliceRList_analyzed_high_confidence.Rdata")
} else {
    load("02_mySpliceRList_analyzed.Rdata")
    mySpliceRList_hc <- mySpliceRList
}

# convert fold change to numeric format
if(!is.null(opt$foldchange)) {
    opt$foldchange <- as.numeric(opt$foldchange)
}

# export GTF files
if(!is.null(opt$generateGTF)) {
    generateGTF(mySpliceRList_hc, filters=c("geneOK", "isoOK", "expressedGenes", "expressedIso"), scoreMethod="local", useProgressBar=T)
    ## uncomment for ying's data (only used in first version)
    #generateGTF(mySpliceRList_hc, filters=c("geneOK"), scoreMethod="local", useProgressBar=T)
}

# retrieve essential features
transcript_features <- as.data.frame(mySpliceRList_hc$transcript_features)
if(ncol(transcript_features)!=78) {
    cat("\nLooks like you are analyzing data from older version of spliceR")
    cat("\nThis script only works for spliceR version 1.8.0")
    cat("\nQuitting.. done\n\n")
    q()
}
exon_features <- unique(as.data.frame(mySpliceRList_hc$exon_features)[,c(5,6)])
mySpliceRList_hc <- merge(transcript_features, exon_features, by.x="spliceR.isoform_id", by.y="spliceR.isoform_id")

## perform independent filtering similar to that performed in DESeq2
if(!is.null(opt$indFilter)) {
    if(file.info(opt$indFilter)$isdir) {
        ncol <- ncol(mySpliceRList_hc)
        ## independent filtering for isoforms
        pdf("independent_filter_threshold.pdf")
        par(mfrow=c(2,2))
        normCounts <- read.table(sprintf("%s/isoforms.count_matrix", opt$indFilter), header=T)
        myIsoforms <- read.table(sprintf("%s/isoform_exp.diff", opt$indFilter), header=T)
        filterstat <- rowMeans(normCounts)
        test <- myIsoforms$p_value
        theta <- seq(0, 1.0, .01)
        alpha=opt$qvalue
        pBH <- filtered_p(filter=filterstat, test=test, theta=theta, method="BH")
        rejBH <- filtered_R(alpha=alpha, filter=filterstat, test=test, theta=theta, method="BH")
        plot(theta, rejBH, type="l", xlab=expression(theta), ylab="Rejections", main=sprintf("BH cutoff = %0.2f", alpha))
        pBH <- as.data.frame(pBH)
        myIsoforms$spliceR.iso_q_value_filter <- pBH[,as.data.frame(which(rejBH==max(rejBH)))[,1][1]]
        mySpliceRList_hc <- merge(mySpliceRList_hc, myIsoforms[,c(1,ncol(myIsoforms))], by.x="spliceR.isoform_id", by.y="test_id")
        mySpliceRList_hc$spliceR.iso_q_value <- mySpliceRList_hc$spliceR.iso_q_value_filter
        mySpliceRList_hc <- mySpliceRList_hc[,c(1:ncol)]

        ## independent filtering for genes
        normCounts <- read.table(sprintf("%s/genes.count_matrix", opt$indFilter), header=T)
        myGenes <- read.table(sprintf("%s/gene_exp.diff", opt$indFilter), header=T)
        filterstat <- rowMeans(normCounts)
        test <- myGenes$p_value
        theta <- seq(0, 1.0, .01)
        alpha=opt$qvalue
        pBH <- filtered_p(filter=filterstat, test=test, theta=theta, method="BH")
        rejBH <- filtered_R(alpha=alpha, filter=filterstat, test=test, theta=theta, method="BH")
        plot(theta, rejBH, type="l", xlab=expression(theta), ylab="Rejections", main=sprintf("BH cutoff = %0.2f", alpha))
        pBH <- as.data.frame(pBH)
        myGenes$spliceR.gene_q_value_filter <- pBH[,as.data.frame(which(rejBH==max(rejBH)))[,1][1]]
        mySpliceRList_hc$spliceR.gene_id_format <- gsub(":.*", "", mySpliceRList_hc$spliceR.gene_id)
        mySpliceRList_hc <- merge(mySpliceRList_hc, myGenes[,c(1,ncol(myGenes))], by.x="spliceR.gene_id_format", by.y="test_id")
        mySpliceRList_hc$spliceR.gene_q_value <- mySpliceRList_hc$spliceR.gene_q_value_filter
        mySpliceRList_hc <- mySpliceRList_hc[,c(2:(ncol+1))]
        dev.off()
    }
}

#if(opt$organism=="mm9") {
#    if(is.null(opt$xlsFormat)) {
#        mySpliceRList_hc$spliceR.gene_id <- sprintf("\"=HYPERLINK(\"\"http://seqbox.local/cgi-bin/hgTracks?clade=mammal&org=Mouse&db=mm9&position=%s:%s-%s\"\",\"\"%s\"\")\"", mySpliceRList_hc$seqnames, mySpliceRList_hc$start, mySpliceRList_hc$end, mySpliceRList_hc$spliceR.gene_id)
#    }
#} else {
#    if(is.null(opt$xlsFormat)) {
#        mySpliceRList_hc$spliceR.gene_id <- sprintf("\"=HYPERLINK(\"\"http://seqbox.local/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg19&position=%s:%s-%s\"\",\"\"%s\"\")\"", mySpliceRList_hc$seqnames, mySpliceRList_hc$start, mySpliceRList_hc$end, mySpliceRList_hc$spliceR.gene_id)
#    }
#}

#mySpliceRList_hc <- mySpliceRList_hc[,c(1,9,11,7,8,2,3,4,5,68,14:length(colnames(mySpliceRList_hc)))] ## all columns excluding 6,10,12,13 (for older version of spliceR)
mySpliceRList_hc <- mySpliceRList_hc[,c(1,9,11,7,8,2,3,4,5,79,6,10,12:length(colnames(mySpliceRList_hc)))]
mySpliceRList_hc <- mySpliceRList_hc[order(-abs(mySpliceRList_hc$spliceR.gene_log2_fold_change), mySpliceRList_hc$spliceR.gene_q_value),]

## add column to data frame describing 5' UTR and 3' UTR length
mySpliceRList_hc$spliceR.utr5Length <- mySpliceRList_hc$spliceR.cdsPosTranscipt
mySpliceRList_hc$spliceR.utr3Length <- mySpliceRList_hc$spliceR.stopPosTranscipt
mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.utr5Length <- (mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.cdsPosTranscipt-1)
mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.utr3Length <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.length-mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.stopPosTranscipt+1

## add column to data frame describing if a transcript is specific to WT or KO
## NOTE: for MIWI data analysis, I considered transcripts having expression <5 and >5 in KO and WT, respectively as specific to WT. This is now changed to <0.5 and >2 for Matilda dataset.
mySpliceRList_hc$spliceR.specificExpression <- "both"
mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1<0.5 & mySpliceRList_hc$spliceR.iso_value_2>opt$minexpr),]$spliceR.specificExpression <- "KO"
mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>opt$minexpr & mySpliceRList_hc$spliceR.iso_value_2<0.5),]$spliceR.specificExpression <- "WT"

## add column to data frame describing if a gene is up-, down- and non-regulated (fold change = KO (sample 2)/WT (sample 1))
mySpliceRList_hc$spliceR.geneRegulation <- "neutral"
if(!is.null(opt$foldchange)) {
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_value_1>opt$minexpr & mySpliceRList_hc$spliceR.gene_value_2>opt$minexpr & mySpliceRList_hc$spliceR.gene_log2_fold_change > opt$foldchange),]$spliceR.geneRegulation <- "up"
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_value_1>opt$minexpr & mySpliceRList_hc$spliceR.gene_value_2>opt$minexpr & mySpliceRList_hc$spliceR.gene_log2_fold_change < -opt$foldchange),]$spliceR.geneRegulation <- "down"
} else if(length(which(mySpliceRList_hc$spliceR.gene_q_value<opt$qvalue)>0)) {
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_q_value<opt$qvalue & mySpliceRList_hc$spliceR.gene_log2_fold_change > 0),]$spliceR.geneRegulation <- "up"
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_q_value<opt$qvalue & mySpliceRList_hc$spliceR.gene_log2_fold_change < 0),]$spliceR.geneRegulation <- "down"
}

## add column to data frame describing if a transcript (isoform) is up-, down- and non-regulated (fold change = KO (sample 2)/WT (sample 1))
mySpliceRList_hc$spliceR.isoRegulation <- "neutral"
if(!is.null(opt$foldchange)) {
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>opt$minexpr & mySpliceRList_hc$spliceR.iso_value_2>opt$minexpr & mySpliceRList_hc$spliceR.iso_log2_fold_change > opt$foldchange),]$spliceR.isoRegulation <- "up"
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>opt$minexpr & mySpliceRList_hc$spliceR.iso_value_2>opt$minexpr & mySpliceRList_hc$spliceR.iso_log2_fold_change < -opt$foldchange),]$spliceR.isoRegulation <- "down"
} else if(length(which(mySpliceRList_hc$spliceR.iso_q_value<opt$qvalue)>0)) {
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_q_value<opt$qvalue & mySpliceRList_hc$spliceR.iso_log2_fold_change > 0),]$spliceR.isoRegulation <- "up"
    mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_q_value<opt$qvalue & mySpliceRList_hc$spliceR.iso_log2_fold_change < 0),]$spliceR.isoRegulation <- "down"
}

# determine transcripts showing switch in their expression pattern between WT and KO
#mySpliceRList_hc_switching <- mySpliceRList_hc[,c(1,6,7,8,2,3,37,38,39)] ## for older version of spliceR
mySpliceRList_hc_switching <- mySpliceRList_hc[,c(1,6,7,8,2,3,29:34,36:40,52,53,54,55,56,57,58,59,60)]
mySpliceRList_hc_switching$spliceR.uid <- paste(mySpliceRList_hc_switching[,5], mySpliceRList_hc_switching[,6], sep="_")
# replaced & with |
mySpliceRList_hc_switching <- mySpliceRList_hc_switching[which((mySpliceRList_hc_switching$spliceR.IF1>opt$minisoexpr | mySpliceRList_hc_switching$spliceR.IF2>opt$minisoexpr) & (mySpliceRList_hc_switching$spliceR.dIF < -opt$minisoexprdiff | mySpliceRList_hc_switching$spliceR.dIF > opt$minisoexprdiff)),]
negative <- mySpliceRList_hc_switching[which(mySpliceRList_hc_switching$spliceR.dIF<0),]
negative <- negative[!duplicated(negative$spliceR.uid),]
positive <- mySpliceRList_hc_switching[which(mySpliceRList_hc_switching$spliceR.dIF>0),]
positive <- positive[!duplicated(positive$spliceR.uid),]
#head(positive[duplicated(positive[order(positive$spliceR.uid, positive$spliceR.dIF),]$spliceR.uid) | duplicated(positive[order(positive$spliceR.uid, positive$spliceR.dIF),]$spliceR.uid, fromLast=TRUE),])
mySpliceRList_hc_switching <- rbind(negative, positive)
mySpliceRList_hc_switching <- mySpliceRList_hc_switching[order(mySpliceRList_hc_switching$spliceR.uid),]
mySpliceRList_hc_switching <- (mySpliceRList_hc_switching[duplicated(mySpliceRList_hc_switching$spliceR.uid) | duplicated(mySpliceRList_hc_switching$spliceR.uid, fromLast = TRUE),])

# determine transcripts showing any of the four splicing events (ESI, MEE, MESI, ISI, A5, A3).
#mySpliceRList_hc_splicing <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.isoRegulation!="neutral"),c(1:3, 6:8, 29:34, 36:40, 52:58, 63)]
#mySpliceRList_hc_splicing <- mySpliceRList_hc[,c(1:3, 6:8, 29:34, 36:40, 52:58, 63)]
# replaced & with |
mySpliceRList_hc_splicing <- mySpliceRList_hc[which((mySpliceRList_hc$spliceR.IF1>opt$minisoexpr | mySpliceRList_hc$spliceR.IF2>opt$minisoexpr) & mySpliceRList_hc$spliceR.isoRegulation!="neutral"),c(1:3, 6:8, 29:34, 36:40, 52:60, 63:75, 80)]

mySpliceRList_hc_splicing$spliced_status <- 0
mySpliceRList_hc_splicing[which(mySpliceRList_hc_splicing$spliceR.ESI>0 | mySpliceRList_hc_splicing$spliceR.MEE>0 | mySpliceRList_hc_splicing$spliceR.MESI>0 | mySpliceRList_hc_splicing$spliceR.ISI>0 | mySpliceRList_hc_splicing$spliceR.A5>0 | mySpliceRList_hc_splicing$spliceR.A3>0),]$spliced_status <- 1
mySpliceRList_hc_splicing <- mySpliceRList_hc_splicing[which(mySpliceRList_hc_splicing$spliced_status==1),]

# differentially expressed genes
#mySpliceRList_hc_gene <- unique(mySpliceRList_hc[,c(2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,69)])
mySpliceRList_hc_gene <- unique(mySpliceRList_hc[,c(2,3,4,5,6,7,8,9,10,28,29,30,31,32,33,34,84)])

if(is.null(opt$xlsFormat)) {
    ## all genes
    write.table(mySpliceRList_hc_gene, file="mySpliceRList_hc_gene.xls", row.names=F, quote=F, sep="\t")
    ## DE genes
    write.table(mySpliceRList_hc_gene[which(mySpliceRList_hc_gene$spliceR.geneRegulation!="neutral"),], file="mySpliceRList_hc_gene_sig.xls", row.names=F, quote=F, sep="\t")
    ## all isoforms
    write.table(mySpliceRList_hc, file="mySpliceRList_hc_iso.xls", row.names=F, quote=F, sep="\t")
    # DE isoforms
    write.table(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.isoRegulation!="neutral"),], file="mySpliceRList_hc_iso_sig.xls", row.names=F, quote=F, sep="\t")
    ## splicing
    write.table(mySpliceRList_hc_splicing, file="mySpliceRList_hc_splicing.xls", row.names=F, quote=F, sep="\t")
    # transcripts showing expression switch
    write.table(mySpliceRList_hc_switching, file="mySpliceRList_hc_switching.xls", row.names=F, quote=F, sep="\t")
    # PTC in transcripts
    write.table(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="TRUE"),], file="mySpliceRList_hc_ptc.xls", quote=F, row.names=F, sep="\t")
    ## WT specific transcripts
    write.table(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT"),], file="mySpliceRList_hc_wt.xls", quote=F, row.names=F, sep="\t")
    ## KO specific transcripts
    write.table(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO"),], file="mySpliceRList_hc_ko.xls", quote=F, row.names=F, sep="\t")
} else {
    ## all genes
    WriteXLS("mySpliceRList_hc_gene", ExcelFileName="mySpliceRList_hc_gene.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## DE genes
    data <- mySpliceRList_hc_gene[which(mySpliceRList_hc_gene$spliceR.geneRegulation!="neutral"),]
    WriteXLS("data", ExcelFileName="mySpliceRList_hc_gene_sig.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## all isoforms
    WriteXLS("mySpliceRList_hc", ExcelFileName="mySpliceRList_hc_iso.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## DE isoforms
    data <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.isoRegulation!="neutral"),]
    WriteXLS("data", ExcelFileName="mySpliceRList_hc_iso_sig.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## splicing
    WriteXLS("mySpliceRList_hc_splicing", ExcelFileName="mySpliceRList_hc_splicing.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## transcripts showing expression switch
    WriteXLS("mySpliceRList_hc_switching", ExcelFileName="mySpliceRList_hc_switching.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## PTC in transcripts
    data <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="TRUE"),]
    WriteXLS("data", ExcelFileName="mySpliceRList_hc_ptc.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## WT specific transcripts
    data <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT"),]
    WriteXLS("data", ExcelFileName="mySpliceRList_hc_wt.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
    ## KO specific transcripts
    data <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO"),]
    WriteXLS("data", ExcelFileName="mySpliceRList_hc_ko.xls", row.names=F, col.names=T, BoldHeaderRow=T, FreezeRow=1)
}

## overlap between gene predictions from cuffdiff and DESeq2
if(!is.null(opt$deseq2SigGeneFile)) {
    if(file.exists(opt$deseq2SigGeneFile)) {
        deseq2_gene <- read.csv(opt$deseq2SigGeneFile, header=T)
        pdf("overlap_gene_deseq2_cuffdiff.pdf")
        venn.plot <- venn.diagram(list(deseq2_gene$gene_short_name, mySpliceRList_hc_gene[which(mySpliceRList_hc_gene$spliceR.geneRegulation!="neutral"),]$spliceR.gene_short_name), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("DESeq2", "Cuffdiff"), main="Differentially expressed genes")
        grid.draw(venn.plot)
        dev.off()
    }
}

## overlap between differentially expressed genes, isoforms, splicing and switching
genes <- gsub(")", "", gsub("\"", "", gsub("^.*,", "", mySpliceRList_hc_gene[which(mySpliceRList_hc_gene$spliceR.geneRegulation!="neutral"),]$spliceR.gene_id)))

isoforms <- unique(gsub(")", "", gsub("\"", "", gsub("^.*,", "", mySpliceRList_hc[which(mySpliceRList_hc$spliceR.isoRegulation!="neutral"),]$spliceR.gene_id))))

switching <- unique(gsub(")", "", gsub("\"", "", gsub("^.*,", "", mySpliceRList_hc_switching$spliceR.gene_id))))

splicing <- unique(gsub(")", "", gsub("\"", "", gsub("^.*,", "", mySpliceRList_hc_splicing$spliceR.gene_id))))

#intersect(intersect(intersect(genes, isoforms), switching), splicing)

pdf("overlap_genes_iso_switching_splicing.pdf")
venn.plot <- venn.diagram(list(genes, isoforms, switching, splicing), NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.2,0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Genes (DE)", "Isoforms (DE)", "Switching (dIF)", "Splicing (ESI, MEE, MESI, ISI, A5, A3)"), main="")
grid.draw(venn.plot)
dev.off()

## overlap between isoform predictions from cuffdiff and DESeq2
if(!is.null(opt$deseq2SigIsoFile)) {
    if(file.exists(opt$deseq2SigIsoFile)) {
        deseq2_iso <- read.table(pipe(sprintf("cut -f 2 -d , %s | sed 's/\"//g' | sed 's/)//g' | grep -v gene_id", opt$deseq2SigIsoFile)))
        pdf("overlap_iso_deseq2_cuffdiff.pdf")
        venn.plot <- venn.diagram(list(deseq2_iso$V1, mySpliceRList_hc[which(mySpliceRList_hc$spliceR.isoRegulation!="neutral"),]$spliceR.isoform_id), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("DESeq2", "Cuffdiff"), main="Differentially expressed isoforms")
        grid.draw(venn.plot)
        dev.off()
    }
}

## overlap between gene and isoform predictions from cuffdiff and deseq2 based on relaxed p-value cut-off of cuffdiff
if(!is.null(opt$deseq2SigGeneFile) && !is.null(opt$deseq2SigIsoFile)) {
    pdf("overlap_pvalue_deseq2_cuffdiff.pdf")
    par(mfrow=c(2,3))
    max_overlap=length(intersect(deseq2_gene$gene_short_name, mySpliceRList_hc$spliceR.gene_short_name))
    j=1;
    overlap<-NULL;
    for(i in seq(0,1,by=0.05)){
        overlap[j]=length(intersect(deseq2_gene$gene_short_name, mySpliceRList_hc_gene[(which(mySpliceRList_hc_gene$spliceR.gene_q_value<i)),]$spliceR.gene_short_name));
        j=j+1;
    }
    barplot(overlap, names.arg=seq(0,1,by=0.05), las=2, ylim=c(0,max_overlap), main="gene", xlab="q-value")
    hist(merge(mySpliceRList_hc_gene, deseq2_gene, by.x="spliceR.gene_short_name", by.y="gene_short_name")$spliceR.gene_p_value, xlab="p-value", main="gene")
    hist(merge(mySpliceRList_hc_gene, deseq2_gene, by.x="spliceR.gene_short_name", by.y="gene_short_name")$spliceR.gene_q_value, xlab="q-value", main="gene")

    max_overlap=length(intersect(deseq2_iso$V1, mySpliceRList_hc$spliceR.isoform_id))
    j=1;
    overlap<-NULL;
    for(i in seq(0,1,by=0.05)){
        overlap[j]=length(intersect(deseq2_iso$V1, mySpliceRList_hc[(which(mySpliceRList_hc$spliceR.iso_q_value<i)),]$spliceR.isoform_id));
        j=j+1;
    }
    barplot(overlap, names.arg=seq(0,1,by=0.05), las=2, ylim=c(0,max_overlap), main="isoform", xlab="q-value")
    hist(merge(mySpliceRList_hc, deseq2_iso, by.x="spliceR.isoform_id", by.y="V1")$spliceR.iso_p_value, xlab="p-value", main="isoform")
    hist(merge(mySpliceRList_hc, deseq2_iso, by.x="spliceR.isoform_id", by.y="V1")$spliceR.iso_q_value, xlab="q-value", main="isoform")
    dev.off()
}

# determine if PTC is significantly enriched in down or up-regulated transcripts
p.value <- fisher.test(matrix(c(length(which(mySpliceRList_hc$spliceR.isoRegulation=="down" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.isoRegulation!="down" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.isoRegulation=="down" & mySpliceRList_hc$spliceR.PTC=="FALSE")), length(which(mySpliceRList_hc$spliceR.isoRegulation!="down" & mySpliceRList_hc$spliceR.PTC=="FALSE"))), nrow=2, byrow=T), alternative = "g")$p.value
cat(sprintf("Significance level at which PTC is enriched in down-regulated transcripts (%d): %02f", length(which(mySpliceRList_hc$spliceR.isoRegulation=="down" & mySpliceRList_hc$spliceR.PTC=="TRUE")), p.value))
cat("\n")

p.value <- fisher.test(matrix(c(length(which(mySpliceRList_hc$spliceR.isoRegulation=="up" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.isoRegulation!="up" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.isoRegulation=="up" & mySpliceRList_hc$spliceR.PTC=="FALSE")), length(which(mySpliceRList_hc$spliceR.isoRegulation!="up" & mySpliceRList_hc$spliceR.PTC=="FALSE"))), nrow=2, byrow=T), alternative = "g")$p.value
cat(sprintf("Significance level at which PTC is enriched in up-regulated transcripts (%d): %02f", length(which(mySpliceRList_hc$spliceR.isoRegulation=="up" & mySpliceRList_hc$spliceR.PTC=="TRUE")), p.value))
cat("\n")

# determine if PTC is significantly enriched in KO or WT specific transcripts
p.value <- fisher.test(matrix(c(length(which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.specificExpression!="KO" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.PTC=="FALSE")), length(which(mySpliceRList_hc$spliceR.specificExpression!="KO" & mySpliceRList_hc$spliceR.PTC=="FALSE"))), nrow=2, byrow=T), alternative = "g")$p.value
cat(sprintf("Significance level at which PTC is enriched in KO specific transcripts (%d): %02f", length(which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.PTC=="TRUE")), p.value))
cat("\n")

p.value <- fisher.test(matrix(c(length(which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.specificExpression!="WT" & mySpliceRList_hc$spliceR.PTC=="TRUE")), length(which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.PTC=="FALSE")), length(which(mySpliceRList_hc$spliceR.specificExpression!="WT" & mySpliceRList_hc$spliceR.PTC=="FALSE"))), nrow=2, byrow=T), alternative = "g")$p.value
cat(sprintf("Significance level at which PTC is enriched in WT specific transcripts (%d): %02f", length(which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.PTC=="TRUE")), p.value))
cat("\n")

## save the session for further use
save.session("mySpliceRList_hc.session")
q()



## old code
## perform independent filtering, if results from DESeq2 are provided (old version)
#if(!is.null(opt$deseq2ResFile)) {
#    if(file.exists(opt$deseq2ResFile)) {
#        deseq2 <- read.table(pipe(sprintf("grep -w NA %s | cut -f 1 -d \",\"", opt$deseq2ResFile)))
#        mySpliceRList_hc$spliceR.gene_id_format <- gsub(":.*", "", mySpliceRList_hc$spliceR.gene_id)
#        mySpliceRList_hc <- mySpliceRList_hc[!mySpliceRList_hc$spliceR.gene_id_format %in% deseq2$V1,c(1:(ncol(mySpliceRList_hc)-1))]
#        mySpliceRList_hc$spliceR.gene_q_value <- p.adjust(mySpliceRList_hc$spliceR.gene_p_value, method="BH")
#        mySpliceRList_hc$spliceR.iso_q_value <- p.adjust(mySpliceRList_hc$spliceR.iso_p_value, method="BH")
#    }
#}

# down or up regulated genes (fold change = KO (sample 2)/WT (sample 1))
#if(!is.null(opt$foldchange)) {
#  mySpliceRList_hc_gene <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_value_1>opt$minexpr & mySpliceRList_hc$spliceR.gene_value_2>opt$minexpr & (mySpliceRList_hc$spliceR.gene_log2_fold_change < -opt$foldchange | mySpliceRList_hc$spliceR.gene_log2_fold_change > opt$foldchange)),]
#} else {
#  mySpliceRList_hc_gene <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.gene_q_value<opt$qvalue & (mySpliceRList_hc$spliceR.gene_log2_fold_change < 0 | mySpliceRList_hc$spliceR.gene_log2_fold_change > 0)),]
#}
# #write.table(sprintf("# min expression: %d; fold change: %d", opt$minexpr, opt$foldchange), file="mySpliceRList_hc_gene.csv", row.names=F, quote=F, col.names=F)
#write.csv(mySpliceRList_hc_gene, file="mySpliceRList_hc_gene.csv", row.names=F, quote=F)
# 
# # down or up regulated isoforms (fold change = KO (sample 2)/WT (sample 1))
# if(!is.null(opt$foldchange)) {
#   mySpliceRList_hc_iso <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>opt$minexpr & mySpliceRList_hc$spliceR.iso_value_2>opt$minexpr & (mySpliceRList_hc$spliceR.iso_log2_fold_change < -opt$foldchange | mySpliceRList_hc$spliceR.iso_log2_fold_change > opt$foldchange)),]
# } else {
#   mySpliceRList_hc_iso <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_q_value<opt$qvalue & (mySpliceRList_hc$spliceR.iso_log2_fold_change < 0 | mySpliceRList_hc$spliceR.iso_log2_fold_change > 0)),]
# }
# #write.table(sprintf("# min expression: %d; fold change: %d", opt$minexpr, opt$foldchange), file="mySpliceRList_hc_iso.csv", row.names=F, quote=F, col.names=F)
# write.csv(mySpliceRList_hc_iso, file="mySpliceRList_hc_iso.csv", row.names=F, quote=F)
