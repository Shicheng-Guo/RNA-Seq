#!/bin/bash

#### usage ####
usage() {
	echo Program: "sra2fastq.sh (convert reads from SRA to FAST[Q/A] format)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: sra2fastq.sh -i <file> [OPTIONS]"
	echo " -i <file>   [reads in SRA format]"
	echo "[OPTIONS]"
	echo " -f          [convert to FASTA. Default is FASTQ]"
    echo " -c          [reads are in color space format (SOLID sequencing)]"
    echo " -p          [reads are from paired-end sequencing]"
	echo " -h          help"
	echo
	exit 0
}

#### parse options ####
while getopts i:fcph ARG; do
	case "$ARG" in
		i) INFILE=$OPTARG;;
		f) TOFASTA=1;;
        c) COLORSPACE=1;;
        p) PAIREDEND=1;;
		h) HELP=1;;
	esac
done

if [ ! $INFILE ]; then
	echo
	usage
fi

ID=`echo $INFILE | sed 's/\..*//g' | sed 's/^.*\///g'`;

## convert SRA to FASTQ
if [ -z "$COLORSPACE" ]; then
    if [ -z "$PAIREDEND" ]; then
        fastq-dump $INFILE
    else
        fastq-dump --split-3 $INFILE
    fi
else
    abi-dump.2.3.5.2 -A $ID $INFILE
    /home/users/sachin/software/bfast-0.7.0a/scripts/solid2fastq -o $ID $ID"_F3.csfasta" $ID"_F3_QV.qual"
    #fastq-dump $INFILE
    #fastqutils csencode $ID.fastq > $ID.fastq.tmp
    #mv $ID.fastq.tmp $ID.fastq
fi
 
## convert FASTQ to FASTA
if [ ! -z "$TOFASTA" ]; then
    fastqutils tofasta $ID.fastq > $ID.fasta
    rm $ID.fastq
fi

exit
