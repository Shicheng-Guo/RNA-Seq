#!/bin/bash
#PBS -l nodes=1:ppn=4

## DEPENDENCIES
READDIR="."
MINLEN=18
MINQUAL=20
FIRSTBASE=1

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
	echo "             [if argument is a fasta file -> fastq-mcf will be executed]"
	echo "             [if argument is adapter sequence -> FastX will be executed]"
    echo "             [if argument is auto -> adapters sequences will be auto-determined using fastqc]"
    echo "             [if argument not provided -> only minimum length, quality and trimming of first few nucleotides will be done]"
	echo " -t          [also output quality filtered reads in FASTA format]"
    echo " -l <int>    [minimum length of the reads (default: 18)]"
    echo " -u <int>    [minimum quality of the reads (default: 20)]"
    echo " -k <int>    [first base to keep (default: 1)]"
    echo " -x <int>    [last base to keep (default: entire read)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:r:a:tl:u:k:x:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		r) READDIR=$OPTARG;;
		a) ADAPTER=$OPTARG;;
		t) FASTA=1;;
        l) MINLEN=$OPTARG;;
        u) MINQUAL=$OPTARG;;
        k) FIRSTBASE=$OPTARG;;
        x) LASTBASE=$OPTARG;;
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
        if [ ! -z "$LASTBASE" ]; then
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -l $LASTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        else
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        fi

        # trim adapters and other artifacts using fastq-mcf
        fastq-mcf -o $READDIR/$ID.clipped.fastq.tmp -l $MINLEN -q $MINQUAL -w 4 -x 10 -t 0 $ADAPTER $READDIR/$ID.clipped.fastq
        mv $READDIR/$ID.clipped.fastq.tmp $READDIR/$ID.clipped.fastq

        # convert fastq to fasta
		if [ ! -z "$FASTA" ]; then
            fastq_to_fasta -Q 33 -n -i $READDIR/$ID.clipped.fastq -o $READDIR/$ID.clipped.fasta 
        fi

		## USE FASTX ##
		#COMMAND_CLIP=`zgrep -v "^>" $ADAPTER | perl -ane 'chomp($_); $_=~s/[^a-zA-Z]+//g; $command.=sprintf("fastx_clipper -a %s -l 18 -Q 33 -n |", $_); END { $command=~s/\|$//g; print $command."\n"; }'`;
		#if [ ! -z "$FASTA" ]; then
		#	zless $FASTQ | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.fastq.clipped
		#	zless $READDIR/$ID.fastq.clipped | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#else
		#	zless $FASTQ | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#fi
	elif [ "$ADAPTER" == "auto" ]; then
        if [ ! -f "$READDIR/$ID"_fastqc"/fastqc_data.txt" ]; then
            qualityCheck.sh -i $FASTQ -q $READDIR
        fi
        zless $READDIR/$ID"_fastqc"/fastqc_data.txt | perl -ane 'if($_=~/\>\>Overrepresented/) { $start=1; } elsif($_=~/END\_MODULE/) { $start=0; } if($start) { print $_; }' | grep -v Overrepresented  | grep -v Sequence | cut -f 1 | sort | uniq | perl -ane '$i++; print ">adapter$i\n$_";' > $ID.ADAPTERS

        # trim first few nucleotides
        if [ ! -z "$LASTBASE" ]; then
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -l $LASTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        else
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        fi

        # trim adapters and other artifacts using fastq-mcf
        fastq-mcf -o $READDIR/$ID.clipped.fastq.tmp -l $MINLEN -q $MINQUAL -w 4 -x 10 -t 0 $ID.ADAPTERS $READDIR/$ID.clipped.fastq
        mv $READDIR/$ID.clipped.fastq.tmp $READDIR/$ID.clipped.fastq

        # convert fastq to fasta
		if [ ! -z "$FASTA" ]; then
            fastq_to_fasta -Q 33 -n -i $READDIR/$ID.clipped.fastq -o $READDIR/$ID.clipped.fasta 
        fi
	else
        # trim first few nucleotides
        if [ ! -z "$LASTBASE" ]; then
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -l $LASTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        else
            zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -m $MINLEN -o $READDIR/$ID.clipped.fastq
        fi

        # trim adapters and other artifacts using FastX
        fastx_clipper -a $ADAPTER -l $MINLEN -Q 33 -i $READDIR/$ID.clipped.fastq | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.clipped.fastq.tmp
        mv $READDIR/$ID.clipped.fastq.tmp $READDIR/$ID.clipped.fastq

        # convert fastq to fasta
		if [ ! -z "$FASTA" ]; then
            fastx_collapser -Q 33 -i $READDIR/$ID.clipped.fastq | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.clipped.fasta
        fi

        # use fastx_clipper with '-n' parameter, if fasta file is empty
		if [ ! -s "$READDIR/$ID.clipped.fastq" ]; then
            # remove beginning of the read until the 'N'
            # trim adapters and other artifacts using FastX
            fastx_clipper -a $ADAPTER -l $MINLEN -Q 33 -i $READDIR/$ID.clipped.fastq -n | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.clipped.fastq.tmp
            mv $READDIR/$ID.clipped.fastq.tmp $READDIR/$ID.clipped.fastq

            # convert fastq to fasta
		    if [ ! -z "$FASTA" ]; then
                zless $READDIR/$ID.clipped.fastq | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>='$MINLEN') { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.clipped.fasta
            fi
		fi
	fi
else
    # trim first few nucleotides
    if [ ! -z "$LASTBASE" ]; then
        zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -l $LASTBASE -m $MINLEN | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 -o $READDIR/$ID.clipped.fastq
    else
        zless $FASTQ | fastx_trimmer -Q 33 -f $FIRSTBASE -m $MINLEN | fastq_quality_trimmer -t $MINQUAL -l $MINLEN -Q 33 | fastx_artifacts_filter -Q 33 -o $READDIR/$ID.clipped.fastq
    fi
fi

## compress output files
gzip $READDIR/$ID.clipped.fastq

if [ ! -z "$FASTA" ]; then
    gzip $READDIR/$ID.clipped.fasta
fi

echo "done"
