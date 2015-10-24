#!/bin/bash
#PBS -l nodes=1:ppn=4

GENOME="mm9"

#### usage ####
usage() {
    echo
	echo Program: "subSampleFasta.sh (retrieve nucleotide sequences corresponding to a transcript id from fasta file)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: subSampleFasta.sh -i <file> -j <string> [OPTIONS]"
	echo "Options:"
    echo " -i <file>   [input FASTA file]"
    echo " -j <file>   [a unique id corresponding to fasta sequence]"
    echo "[OPTIONS]"
    echo " -g <genome> [reference genome (default: mm9)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:j:g:h ARG; do
	case "$ARG" in
        i) FASTAFILE=$OPTARG;;
        j) FASTAID=$OPTARG;;
        g) GENOME=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$FASTAFILE" -o -z "$FASTAID" -o "$HELP" ]; then
	usage
fi

###################
#helperfunction
function wait_for_jobs_to_finish {
    for job in `jobs -p`
    do
        echo $job
        wait $job
    done
    echo $1
}
###############

if [ "$GENOME" == "mm9" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/mouse.mm9.fa"
elif [ "$GENOME" == "hg19" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/human.hg19.fa"
else
    echo "Presently the program only support analysis for mm9 or hg19"
echo
    usage
fi

if [ ! -f "$FASTAFILE.fai" ]; then
    samtools faidx $FASTAFILE
fi

samtools faidx $FASTAFILE $FASTAID

