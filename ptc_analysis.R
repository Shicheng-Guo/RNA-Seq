wd<- c("/home/pundhir/project/rna-seq-analysis/matilda/result_PreMegE/output/5_spliceR_analysis", "/home/pundhir/project/rna-seq-analysis/matilda/result_CFUE/output/5_spliceR_analysis")
wd_des <- c("rspd", "spc")

suppressPackageStartupMessages(library("spliceR"))

## initialize data frame to store all results
all <- data.frame()

for (i in seq(1, length(wd), 1)) {
  setwd(wd[i])
  load("03_mySpliceRList_analyzed_high_confidence.Rdata")

  transcript_features <- as.data.frame(mySpliceRList_hc$transcript_features)
  exon_features <- unique(as.data.frame(mySpliceRList_hc$exon_features)[,c(5,6)])
  mySpliceRList_hc <- merge(transcript_features, exon_features, by.x="spliceR.isoform_id", by.y="spliceR.isoform_id")

  ## add column to data frame describing 5' UTR and 3' UTR length
  mySpliceRList_hc$spliceR.utr5Length <- mySpliceRList_hc$spliceR.cdsPosTranscipt
  mySpliceRList_hc$spliceR.utr3Length <- mySpliceRList_hc$spliceR.stopPosTranscipt
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.utr5Length <- (mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.cdsPosTranscipt-1)
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.utr3Length <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.length-mySpliceRList_hc[which(mySpliceRList_hc$spliceR.cdsPosTranscipt!=-1),]$spliceR.stopPosTranscipt+1

  ## add column to data frame describing if a transcript is specific to WT or KO
  mySpliceRList_hc$spliceR.specificExpression <- "both"
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1<5 & mySpliceRList_hc$spliceR.iso_value_2>5),]$spliceR.specificExpression <- "KO"
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>5 & mySpliceRList_hc$spliceR.iso_value_2<5),]$spliceR.specificExpression <- "WT"

  ## add column to data frame describing if a transcript is up-, down- and non-regulated
  mySpliceRList_hc$spliceR.regulation <- "neutral"
  #mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>5 & mySpliceRList_hc$spliceR.iso_value_2>5 & (mySpliceRList_hc$spliceR.iso_log2_fold_change < 2 & mySpliceRList_hc$spliceR.iso_log2_fold_change > -2)),]$spliceR.regulation <- "neutral"
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>5 & mySpliceRList_hc$spliceR.iso_value_2>5 & mySpliceRList_hc$spliceR.iso_log2_fold_change > 2),]$spliceR.regulation <- "up"
  mySpliceRList_hc[which(mySpliceRList_hc$spliceR.iso_value_1>5 & mySpliceRList_hc$spliceR.iso_value_2>5 & mySpliceRList_hc$spliceR.iso_log2_fold_change < -2),]$spliceR.regulation <- "down"

  ## make a data frame for all three sample analysis results (rspd, spc and GSE32183)
  mySpliceRList_hc$spliceR.sample <- wd_des[i]
  all <- rbind(all, mySpliceRList_hc)
}

mySpliceRList_hc <- all[which(all$spliceR.sample=="rspd"),]; setwd(wd[1]);
mySpliceRList_hc <- all[which(all$spliceR.sample=="spc"),]; setwd(wd[2]);


## 5' UTR length distribution for isoforms specific to KO and WT
pdf("rspd_5utr_length.pdf")
library(ggplot2)
library(plyr)
ggplot(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression!="both" & mySpliceRList_hc$spliceR.utr5Length!=-1),], aes(x=log2(spliceR.utr5Length+1), fill=spliceR.specificExpression)) + geom_density(alpha=.3) + xlab("5' UTR length (log2)")
dev.off()

wt <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length
ko <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length
round(wilcox.test(ko, wt, alternative="g")$p.value, digits=6)

# median 5' UTR length in WT and KO transcripts
median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length)
median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length)

## 3' UTR length distribution for isoforms specific to KO and WT
pdf("rspd_3utr_length.pdf")
library(ggplot2)
library(plyr)
ggplot(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression!="both" & mySpliceRList_hc$spliceR.utr5Length!=-1),], aes(x=log2(spliceR.utr3Length+1), fill=spliceR.specificExpression)) + geom_density(alpha=.3) + xlab("3' UTR length (log2)")
dev.off()

wt <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length
ko <- mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length
round(wilcox.test(ko, wt, alternative="g")$p.value, digits=100)

# median 3' UTR length in WT and KO transcripts
median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="WT" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length)
median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.specificExpression=="KO" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length)


## plot 5' UTR length with respect to up-, neutral- and down-regulated genes
pdf("rspd_5utr_length_expression.pdf")
mySpliceRList_hc_up <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr5Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
mySpliceRList_hc_neutral <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr5Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
mySpliceRList_hc_down <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr5Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
barplot(rbind(mySpliceRList_hc_down, mySpliceRList_hc_neutral, mySpliceRList_hc_up), beside=1, col=c("blue", "green", "red"), names.arg=(seq(0,15,by=0.5)[2:31]), las=2, xlab="5' UTR length (log2)", ylab="density", cex=1.2, cex.lab=1.2, cex.axis=1.2, cex.names=1.2)
legend("topleft", legend=c(paste("down n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down")), sep=""), paste("neutral n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral")), sep=""), paste("up n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up")), sep="")), fill=c("blue", "green", "red"), cex=1, border="white", bty="n")
dev.off()

## plot 3' UTR length with respect to up-, neutral- and down-regulated genes
pdf("rspd_3utr_length_expression.pdf")
mySpliceRList_hc_up <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr3Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
mySpliceRList_hc_neutral <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr3Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
mySpliceRList_hc_down <- hist(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr3Length+1), plot=FALSE, breaks=seq(0,15,by=0.5))$density/2
barplot(rbind(mySpliceRList_hc_down, mySpliceRList_hc_neutral, mySpliceRList_hc_up), beside=1, col=c("blue", "green", "red"), names.arg=(seq(0,15,by=0.5)[2:31]), las=2, xlab="3' UTR length (log2)", ylab="density", cex=1.2, cex.lab=1.2, cex.axis=1.2, cex.names=1.2)
legend("topleft", legend=c(paste("down n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down")), sep=""), paste("neutral n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral")), sep=""), paste("up n=", length(which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up")), sep="")), fill=c("blue", "green", "red"), cex=1, border="white", bty="n")
dev.off()

# median 5' UTR length in up, down and neutral regulated transcripts
pdf("rspd_5utr_length_expression_pVal.pdf")
median <- c(median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr5Length+1)), median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length+1)), median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length+1)))
mp <- barplot(median, axes=FALSE, axisnames=FALSE, ylim=c(0,18), col=c("red", "green", "blue"), ylab="median 5' UTR length (log2)", cex=1.2, cex.axis=1.2, cex.lab=1.2)
axis(1, labels=c("up", "neutral", "down"), at=mp, cex=1.2, cex.axis=1.2, cex.lab=1.2)
axis(2, at=seq(0,18,by=2), cex.axis=1.2, cex.lab=1.2)
box()
sd <- c(sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr5Length+1)), sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length+1)), sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length+1)))
text(mp[1], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="up" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
text(mp[2], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
text(mp[3], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
text(mp[1], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="up" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
text(mp[2], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
text(mp[3], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr5Length))), cex=1.2)
segments(mp, median-sd, mp, median+sd, lwd=2)
segments(mp-0.1, median - sd, mp + 0.1, median - sd, lwd=2)
segments(mp-0.1, median + sd, mp + 0.1, median + sd, lwd=2)
segments(mp[1], 12, mp[2], 12, lwd=2); text((mp[1]+mp[2])/2, 12+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr5Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr5Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
segments(mp[2], 14, mp[3], 14, lwd=2); text((mp[2]+mp[3])/2, 14+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr5Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr5Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
segments(mp[1], 16, mp[3], 16, lwd=2); text((mp[1]+mp[3])/2, 16+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr5Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr5Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
dev.off()

# median 3' UTR length in up, down and neutral regulated transcripts
pdf("rspd_3utr_length_expression_pVal.pdf")
median <- c(median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr3Length+1)), median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length+1)), median(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length+1)))
mp <- barplot(median, axes=FALSE, axisnames=FALSE, ylim=c(0,18), col=c("red", "green", "blue"), ylab="median 3' UTR length (log2)", cex=1.2, cex.axis=1.2, cex.lab=1.2)
axis(1, labels=c("up", "neutral", "down"), at=mp, cex=1.2, cex.axis=1.2, cex.lab=1.2)
axis(2, at=seq(0,18,by=2), cex.axis=1.2, cex.lab=1.2)
box()
sd <- c(sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr3Length+1)), sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length+1)), sd(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length+1)))
text(mp[1], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="up" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
text(mp[2], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
text(mp[3], 0.8, sprintf("median= %d", round(median(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
text(mp[1], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="up" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
text(mp[2], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="neutral" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
text(mp[3], 2.0, sprintf("mean= %d", round(mean(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.regulation=="down" & mySpliceRList_hc$spliceR.utr5Length!=-1),]$spliceR.utr3Length))), cex=1.2)
segments(mp, median-sd, mp, median+sd, lwd=2)
segments(mp-0.1, median - sd, mp + 0.1, median - sd, lwd=2)
segments(mp-0.1, median + sd, mp + 0.1, median + sd, lwd=2)
segments(mp[1], 12, mp[2], 12, lwd=2); text((mp[1]+mp[2])/2, 12+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr3Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr3Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
segments(mp[2], 14, mp[3], 14, lwd=2); text((mp[2]+mp[3])/2, 14+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="neutral"),]$spliceR.utr3Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr3Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
segments(mp[1], 16, mp[3], 16, lwd=2); text((mp[1]+mp[3])/2, 16+0.6, labels=sprintf("p-value= %s", formatC(wilcox.test(log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="up"),]$spliceR.utr3Length+1), log2(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.utr5Length!=-1 & mySpliceRList_hc$spliceR.regulation=="down"),]$spliceR.utr3Length+1), alternative="g")$p.value, digits=2)), cex=1.2)
dev.off()

## PTC frequency among up-regulated transcripts
length(which(mySpliceRList_hc$spliceR.PTC=="TRUE" & mySpliceRList_hc$spliceR.regulation=="down"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="TRUE"),]$spliceR.PTC)*100
length(which(mySpliceRList_hc$spliceR.PTC=="TRUE" & mySpliceRList_hc$spliceR.regulation=="neutral"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="TRUE"),]$spliceR.PTC)*100
length(which(mySpliceRList_hc$spliceR.PTC=="TRUE" & mySpliceRList_hc$spliceR.regulation=="up"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="TRUE"),]$spliceR.PTC)*100

length(which(mySpliceRList_hc$spliceR.PTC=="FALSE" & mySpliceRList_hc$spliceR.regulation=="down"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="FALSE"),]$spliceR.PTC)*100
length(which(mySpliceRList_hc$spliceR.PTC=="FALSE" & mySpliceRList_hc$spliceR.regulation=="neutral"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="FALSE"),]$spliceR.PTC)*100
length(which(mySpliceRList_hc$spliceR.PTC=="FALSE" & mySpliceRList_hc$spliceR.regulation=="up"))/length(mySpliceRList_hc[which(mySpliceRList_hc$spliceR.PTC=="FALSE"),]$spliceR.PTC)*100


length(mySpliceRList_hc[(which(mySpliceRList_hc$spliceR.specificExpression=="KO" & !is.na(mySpliceRList_hc$spliceR.stopPosGenomic))),]$spliceR.nearest_ref_id)/length(mySpliceRList_hc[(which(mySpliceRList_hc$spliceR.specificExpression=="KO")),]$spliceR.nearest_ref_id)
length(mySpliceRList_hc[(which(mySpliceRList_hc$spliceR.specificExpression=="WT" & !is.na(mySpliceRList_hc$spliceR.stopPosGenomic))),]$spliceR.nearest_ref_id)/length(mySpliceRList_hc[(which(mySpliceRList_hc$spliceR.specificExpression=="WT")),]$spliceR.nearest_ref_id)



########################################################################################################
## determine common and unique KO specific transcripts expressed in rspd, spc and GEO32183 samples
rspd <- all[which(all$spliceR.sample=="rspd" & all$spliceR.specificExpression=="KO" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id
spc <- all[which(all$spliceR.sample=="spc" & all$spliceR.specificExpression=="KO" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id
geo32183 <- all[which(all$spliceR.sample=="geo32183" & all$spliceR.specificExpression=="KO" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id

temp <- venn.diagram(x=list(A=rspd, B=spc, C=geo32183),
             category.names=c("RSPD_KO", "SPC_KO", "GEO32183_KO"),
             filename=NULL,
             scaled=TRUE,
             fill = c("red", "blue", "green"),
             cex=2.5
             )
pdf("KO.pdf")
grid.draw(temp)
dev.off()

## determine common and unique WT specific transcripts expressed in rspd, spc and GEO32183 samples
rspd <- all[which(all$spliceR.sample=="rspd" & all$spliceR.specificExpression=="WT" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id
rspd <- rspd[-278]
spc <- all[which(all$spliceR.sample=="spc" & all$spliceR.specificExpression=="WT" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id
geo32183 <- all[which(all$spliceR.sample=="geo32183" & all$spliceR.specificExpression=="WT" & all$spliceR.utr5Length!=-1),]$spliceR.nearest_ref_id

temp <- venn.diagram(x=list(A=rspd, B=spc, C=geo32183),
                     category.names=c("RSPD_WT", "SPC_WT", "GEO32183_WT"),
                     filename=NULL,
                     scaled=TRUE,
                     fill = c("red", "blue", "green"),
                     cex=2.5
)
pdf("WT.pdf")
grid.draw(temp)
dev.off()

