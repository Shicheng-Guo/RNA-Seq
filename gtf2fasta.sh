#!/bin/bash
#PBS -l nodes=1:ppn=4

GENOME="mm9"

#### usage ####
usage() {
    echo
	echo Program: "gtf2fasta.sh (retrieve nucleotide sequences corresponding to transcripts in gtf format)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: gtf2fasta.sh -i <file> [OPTIONS]"
	echo "Options:"
    echo " -i <file>   [input GTF file]"
    echo " -o <file>   [output fasta file]"
    echo "[OPTIONS]"
    echo " -g <genome> [reference genome (default: mm9)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:o:g:h ARG; do
	case "$ARG" in
		i) GTFFILE=$OPTARG;;
        o) FASTAFILE=$OPTARG;;
        g) GENOME=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$GTFFILE" -o -z "$FASTAFILE" -o "$HELP" ]; then
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

echo -n "Populating files based on input genome ($GENOME)... "
if [ "$GENOME" == "mm9" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/mouse.mm9.fa"
elif [ "$GENOME" == "hg19" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/human.hg19.fa"
else
    echo "Presently the program only support analysis for mm9 or hg19"
echo
    usage
fi
echo done

echo -n "format GTF file... "
ID=$(echo $GTFFILE | sed 's/\.GTF/_mod/ig');

grep -v "track" $GTFFILE | perl -an -F'/\t/' -e '@t=split(/\;/, $F[8]); $line=(); foreach(@t) { @u=split(/\s+/,$_); $des=(); foreach(@u[0..scalar(@u)-2]) { if($_=~/^$/) { next; } $des.=$_."_"; } $des=~s/\_$//g; $line.="$des \"$u[scalar(@u)-1]\"\; "; } print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$line\n";' | perl -an -F'/\t/' -e '@t=split(/\;/,$F[8]); $id="NA"; foreach(@t) { if($_=~/isoform/i || $_=~/transcript/i) { $_=~s/^.*\s+//g; $_=~s/\"//g;  $id=$_;} } print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$id\n";' > $ID.GTF
echo "done"

echo -n "retrieve transcript sequence in fasta format... "
gffread -w $FASTAFILE -g $GENOME_FILE $ID.GTF
rm $ID.GTF
echo "done"

