#!/bin/bash
#PBS -l nodes=1:ppn=4

QUALITYDIR="."

#### usage ####
usage() {
	echo Program: "qualityCheck.sh (perform quality check on input fastq files)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: qualityCheck.sh -i <file> [OPTIONS]"
	echo " -i <file>   [input fastq file with single end reads]"
    echo "[OPTIONS]"
	echo " -q <dir>    [directory to store read quality report (default: .)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:q:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		q) QUALITYDIR=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$FASTQ" -o "$HELP" ]; then
	usage
fi

## check if output directory exists
echo -n "Create appropriate directory structure... "
if [ ! -d "$QUALITYDIR" ]; then
    mkdir -p $QUALITYDIR
fi
echo "done"

## retrieve file name
ID=`echo $FASTQ | perl -ane '$_=~s/^.+\///g; $_=~s/\..+$//g; print $_;'`;

echo -n "Compute quality for $ID... "

eval fastqc \"$FASTQ\" -o $QUALITYDIR

echo "done"
