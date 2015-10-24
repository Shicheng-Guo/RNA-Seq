#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-n", "--genome"), default="hg19", help="reference genome [default %default]"),
	make_option(c("-c", "--coor"), metavar="coordinate", help="coordinate (chr:start-end)"),
	make_option(c("-f", "--file"), metavar="file", help="coordinate file (chr:start-end)"),
	make_option(c("-t", "--trackfile"), help="custom track file"),
	make_option(c("-o", "--output"), help="output pdf file")
)

opt <- parse_args(OptionParser(option_list=option_list))

####### lib ############
## CREDITS
## source: http://biostar.stackexchange.com/questions/6149/batch-viewing-of-ucsc-browser-graphic/6155#6155
## source: http://onetipperday.blogspot.dk/2012/07/get-ucsc-images-for-list-of-regions-in.html
## merge script
mergePDF <- function(output="merged.pdf", sourcefiles=c("source1.pdf","source2.pdf","source3.pdf"))
{
	## create the command string and call the command using system()
	## merging command from http://hints.macworld.com/article.php?story=2003083122212228
	command=paste("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite",paste("-sOutputFile",output, sep="="), paste(sourcefiles, collapse=" "),sep=" ")
	try(system(command))
}

## Here is an R script wrote by Aaron Statham which saves UCSC to pdfs -
## you can choose which genome and tracks to display by altering the 'url' parameter. 'trackfile' is the url of a file describing the custom tracks (beds/bigwigs) to display
screenshotUCSC <- function(url, trackfile, chr, start, end, filename) {
	oldpen <- options("scipen")
	options(scipen=100)
	temp <- readLines(paste(url, "&hgt.customText=", trackfile, "&position=",chr,":",start,"-",end, sep=""))
	pdfurl <- paste("http://genome.ucsc.edu/trash",gsub(".*trash","",gsub(".pdf.*","",temp[grep(".pdf", temp, fixed=TRUE)][1])), ".pdf", sep="")
	options(scipen=oldpen)
	download.file(pdfurl, filename, mode="wb", quiet=TRUE)
}

## define the default track
## controlling individual tracks
theURL="http://genome.ucsc.edu/cgi-bin/hgTracks?db="
theURL=paste(theURL, opt$trackfile,sep="")
theURL=paste(theURL,"&wgRna=hide&cpgIslandExt=pack&ensGene=hide&mrna=hide&intronEst=hide&mgcGenes=hide&hgt.psOutput=on&cons44way=hide&snp130=hide&snpArray=hide&wgEncodeReg=hide&pix=1000&refGene=pack&knownGene=hide&rmsk=hide",sep="")
## controling via session
#theURL="http://genome-preview.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sterding&hgS_otherUserSessionName=<SESSION_NAME>&hgt.psOutput=on&pix=1000"

## read coordinates
toPlot <- NULL
chr <- NULL
start <- NULL
end <- NULL
if ( !is.null(opt$coor) ) {
	toPlot = opt$coor
	chr <- as.vector(as.data.frame(strsplit(opt$coor, "[:-]+"))[1,1])
	start <- as.numeric(as.vector(as.data.frame(strsplit(opt$coor, "[:-]+"))[2,1]))-1000
	end <- as.numeric(as.vector(as.data.frame(strsplit(opt$coor, "[:-]+"))[3,1]))+999
	paste("region_1_", chr, "_", start, "_", end, ".pdf", sep="")
	screenshotUCSC(theURL, "", chr, start, end, paste("region_1_", chr, "_", start, "_", end, ".pdf", sep=""))
} else if ( !is.null(opt$file) ) {
	toPlot=read.table(opt$file)
	## anti-robot version 
	## UCSC Policy: Program-driven use of this software is limited to a maximum of one hit every 15 seconds and no more than 5,000 hits per day.
	for(i in 1:nrow(toPlot)){
		chr <- as.character(toPlot[i,1])
		start <- toPlot[i,2]-1000
		end <- toPlot[i,3]+999
		a<- paste("region_", i, "_", chr, "_", start, "_", end, ".pdf", sep="")
		screenshotUCSC(theURL, "", chr, start, end, paste("region_", i, "_", chr, "_", start, "_", end, ".pdf", sep=""))
		Sys.sleep(5) 
	}
}

mergePDF("piRNAs_ucsc_screenshot.pdf", list.files(pattern="region_.*.pdf"))
try(system("rm region_*.pdf"))

q()
## parallel version
#library(multicore)
#mclapply(1:nrow(toPlot), function(i) screenshotUCSC(theURL, "", as.character(toPlot$chr[i]), toPlot$start[i]-2000, toPlot$end[i]+1999, paste("region_", i, "_", toPlot$name[i],".pdf", sep="")), mc.cores=10)
#"URL_of_your_custom_track", toPlot$space[i], toPlot$start[i]-3000, toPlot$start[i]+2999, paste("Figures/Shots/Low_To_High_", i, ".pdf", sep="")), mc.cores=10)

