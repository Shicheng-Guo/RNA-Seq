#!/bin/bash
#PBS -l nodes=1:ppn=4

## DEPENDENCIES
READDIR="."
MINLEN=18
MINQUAL=20
TRIMLEN=0

#### usage ####
usage() {
	echo Program: "trimAdapters.sh (trim adapter sequences from fastq files)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: trimAdapters.sh -i <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input file with single end reads]"
	echo "             [FAST(Q|A) format - raw file before quality filter]"
	echo "[OPTIONS]"
	echo " -r <dir>    [directory to store quality filtered reads (default: .)]"
	echo " -a <file>   [adapter sequences]"
	echo "             [if argument is a file -> fastq-mcf will be executed]"
	echo "             [if argument is adapter sequence -> FastX will be executed]"
    echo "             [if argument not provided -> only minimum length, quality and trimming of first few nucleotides will be done]"
	echo " -t          [also output quality filtered reads in FASTA format]"
    echo " -l <int>    [minimum length of the reads (default: 18)]"
    echo " -u <int>    [minimum quality of the reads (default: 20)]"
    echo " -k <int>    [trim first few nucleotides (default: 0)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:r:a:tl:u:k:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		r) READDIR=$OPTARG;;
		a) ADAPTER=$OPTARG;;
		t) FASTA=1;;
        l) MINLEN=$OPTARG;;
        u) MINQUAL=$OPTARG;;
        k) TRIMLEN=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$FASTQ" -o "$HELP" ]; then
	usage
fi

## create appropriate directory structure
echo -n "Create appropriate directory structure... "
if [ ! -d "$READDIR" ]; then
    mkdir $READDIR
fi
echo "done"

## retrieve file name
ID=`echo $FASTQ | perl -ane '$_=~s/^.+\///g; $_=~s/\..+$//g; print $_;'`;

echo "perform quality filter (trim adapters etc) for $ID... "

## trim adapters
if [ ! -z "${ADAPTER}" -a ! -z "${READDIR}" ]; then
	if [ -f "$ADAPTER" ]; then
        # trim first few nucleotides
        TRIMLEN=$(($TRIMLEN + 1))
        fastx_trimmer -Q 33 -f $TRIMLEN -i $FASTQ -m $MINLEN -o $READDIR/clipped_$ID.fastq

        # trim adapters and other artifacts using fastq-mcf
        fastq-mcf -o $READDIR/clipped_$ID.fastq.tmp -l $MINLEN -q $MINQUAL -w 4 -x 10 -t 0 $ADAPTER $READDIR/clipped_$ID.fastq
        mv $READDIR/clipped_$ID.fastq.tmp $READDIR/clipped_$ID.fastq

        # convert fastq to fasta
		if [ ! -z "$FASTA" ]; then
            fastq_to_fasta -Q 33 -n -i $READDIR/clipped_$ID.fastq -o $READDIR/clipped_$ID.fasta 
        fi

		## USE FASTX ##
		#COMMAND_CLIP=`zgrep -v "^>" $ADAPTER | perl -ane 'chomp($_); $_=~s/[^a-zA-Z]+//g; $command.=sprintf("fastx_clipper -a %s -l 18 -Q 33 -n |", $_); END { $command=~s/\|$//g; print $command."\n"; }'`;
		#if [ ! -z "$FASTA" ]; then
		#	zless $FASTQ | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.fastq.clipped
		#	zless $READDIR/$ID.fastq.clipped | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#else
		#	zless $FASTQ | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#fi
	else
        # trim first few nucleotides
        TRIMLEN=$(($TRIMLEN + 1))
        fastx_trimmer -Q 33 -f $TRIMLEN -i $FASTQ -m $MINLEN -o $READDIR/clipped_$ID.fastq

        # trim adapters and other artifacts using FastX
        fastx_clipper -a $ADAPTER -l $MINLEN -Q 33 -i $READDIR/clipped_$ID.fastq | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/clipped_$ID.fastq.tmp
        mv $READDIR/clipped_$ID.fastq.tmp $READDIR/clipped_$ID.fastq

        # convert fastq to fasta
		if [ ! -z "$FASTA" ]; then
            fastx_collapser -Q 33 -i $READDIR/clipped_$ID.fastq | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/clipped_$ID.fasta
        fi

        # use fastx_clipper with '-n' parameter, if fasta file is empty
		if [ ! -s "$READDIR/clipped_$ID.fastq" ]; then
            # remove beginning of the read until the 'N'
            # trim adapters and other artifacts using FastX
            fastx_clipper -a $ADAPTER -l $MINLEN -Q 33 -i $READDIR/clipped_$ID.fastq -n | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/clipped_$ID.fastq.tmp
            mv $READDIR/clipped_$ID.fastq.tmp $READDIR/clipped_$ID.fastq

            # convert fastq to fasta
		    if [ ! -z "$FASTA" ]; then
                zless $READDIR/clipped_$ID.fastq | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>='$MINLEN') { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/clipped_$ID.fasta
            fi
		fi
	fi
else
    # trim first few nucleotides
    TRIMLEN=$(($TRIMLEN + 1))
    fastx_trimmer -Q 33 -f $TRIMLEN -i $FASTQ -m $MINLEN | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 -o $READDIR/clipped_$ID.fastq
fi

echo "done"
