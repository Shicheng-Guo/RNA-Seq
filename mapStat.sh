#!/bin/bash
#PBS -l nodes=1:ppn=4

#### usage ####
usage() {
	echo Program: "mapStat (tabulate mapping statistics)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: mapStat.sh -s <file> -q <file> -f <file>"
	echo "Options:"
	echo " -s <file>   [input mapStat file from mapSEReads.sh]"
	echo "[optional]"
	echo " -q <file>   [input FASTQ (raw reads) file]"
	echo " -f <file>   [input FASTQ (clipped read) file]"
    echo " -m <file>   [input BAM file containing mapped reads]"
    echo "             [will be used to output unmapped reads in .unmapped file]"
	echo " -h          [help]"
	echo
	exit 0
}

MAPPING_FREQUENCY=1

#### parse options ####
while getopts s:q:f:m:h ARG; do
	case "$ARG" in
		s) MAPSTATFILE=$OPTARG;;
		q) RAWFASTQFILE=$OPTARG;;
		f) CLIPPEDFASTQFILE=$OPTARG;;
        m) BAMFILE=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$MAPSTATFILE" -o "$HELP" ]; then
	usage
fi

## count total number of raw reads
RAW_READS_COUNT="NA"
if [ ! -z "$RAWFASTQFILE" -a -f "$RAWFASTQFILE" ]; then
    RAW_READS_COUNT=$(zless $RAWFASTQFILE | grep "@" | wc -l)
fi

## count total number of reads left after quality check
CLIPPED_READS_COUNT="NA"
PER="NA"
if [ ! -z "$CLIPPED_READS_COUNT" -a -f "$CLIPPEDFASTQFILE" ]; then
    CLIPPED_READS_COUNT=$(zless $CLIPPEDFASTQFILE | grep "@" | wc -l)
    if [ ! -z "$RAWFASTQFILE" -a -f "$RAWFASTQFILE" ]; then
        PER=$(echo $CLIPPED_READS_COUNT | perl -ane 'chomp($_); $per=($_*100)/'$RAW_READS_COUNT'; printf("%0.2f%", $per);')
    fi
fi

## tabulate mapping statistics
echo -e "fastq_file\t#reads (raw)\t#reads (after qualityCheck)\t#reads (for mapping)\t# reads (unpaired)\t#reads (unmapped)\t#reads (aligned 1 time)\t#reads (aligned >1 time)\talignment rate"
zless $MAPSTATFILE | perl -ane 'BEGIN { print "'$RAWFASTQFILE'\t'$RAW_READS_COUNT'\t'$CLIPPED_READS_COUNT' ('$PER')"; } if($_=~/^[0-9\s]+/) { $_=~s/\;.*//g; chomp($_); $_=~s/\s+[a-zA-Z]+.*//g; print "\t$_"; } END { print "\n"; }'

if [ ! -z "$BAMFILE" ]; then
    ID=`echo $BAMFILE | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
    samtools view -f 4 $BAMFILE  | cut -f 10 | sort | uniq -c | sort -nr > $ID.unmapped
fi

exit
