#!/bin/bash
#PBS -l nodes=1:ppn=4

## DEPENDENCIES
trimSE="java -classpath /home/users/sachin/software/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticSE"
READDIR="."
UCSCDIR="."
MINLEN=18
MINQUAL=20
TRIMLEN=0
PROCESSORS=1
FASTAFILE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/TopHatTranscriptomeIndex/genes_without_mt"
GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/Bowtie2IndexWithAbundance"
CHRSIZE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/ChromInfoRef.txt"

#### usage ####
usage() {
	echo Program: "mapSEReads.sh (perform quqlity check, filter and mapping of \"single-end\" reads)"
	echo Author: RTH, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: sachin@rth.dk
	echo "Usage: mapSEReads.sh -i <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input file with single end reads]"
	echo "             [ FAST(Q|A) format - raw file before quality filter]"
	echo "[Quality check]"
	echo " -q <dir>    [directory to store read quality report]"
	echo "[Trim adapters]"
	echo " -r <dir>    [directory to store quality filtered reads (FASTQ format)]"
	echo " -a <file>   [adapter sequences]"
	echo "             [if argument is a file -> fastq-mcf will be executed]"
	echo "             [if argument is adapter sequence -> FastX will be executed]"
	echo " -t          [also output quality filtered reads in FASTA format]"
    echo " -l <int>    [minimum length of the reads (default: 18)]"
    echo " -u <int>    [minimum quality of the reads (default: 20)]"
    echo " -k <int>    [trim first few nucleotides (default: 0)]"
	echo "[Map reads]"
	echo " -m <dir>    [directory to store mapped reads (SAM format)]"
	echo " -r <dir>    [directory having filtered reads (FASTA format)]"
    echo " -s <dir>    [directory to keep files related to UCSC genome browser (default: current)]"
	echo " -f <file>   [genome file]"
    echo "             [segemehl: file in fasta format]"
    echo "             [tophat:   directory having transcriptome index]"
    echo "             [default: ~/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/TopHatTranscriptomeIndex/genes_without_mt]"
	echo " -g <file>   [indexed genome files]"
    echo "             [segemehl: file in .idx extension]"
    echo "             [tophat:   directory having genome index]"
    echo "             [default: ~/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/Bowtie2IndexWithAbundance]" 
    echo " -c <file>   [file having chromosome's size. required for genomecov]"
    echo "             [default: ~/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/ChromInfoRef.txt]"
    echo " -p <int>    [number of processors (default: 1)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:q:r:a:tl:u:k:m:s:f:g:c:p:h ARG; do
	case "$ARG" in
		i) READ1=$OPTARG;;
		q) QUALITYDIR=$OPTARG;;
		r) READDIR=$OPTARG;;
		a) ADAPTER=$OPTARG;;
		t) FASTA=1;;
        l) MINLEN=$OPTARG;;
        u) MINQUAL=$OPTARG;;
        k) TRIMLEN=$OPTARG;;
		m) MAPDIR=$OPTARG;;
        s) UCSCDIR=$OPTARG;;
		f) FASTAFILE=$OPTARG;;
		g) GENOMEINDEX=$OPTARG;;
        c) CHRSIZE=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$READ1" -o "$HELP" ]; then
	usage
fi

## retrieve file name
ID=`echo $READ1 | perl -ane '$_=~s/^.+\///g; $_=~s/\..+$//g; print $_;'`;

echo "Computing for $ID"

## quality check
if [ ! -z "${QUALITYDIR}" ]; then
	eval fastqc \"$READ1\" -o $QUALITYDIR
fi

## trim adapters
if [ ! -z "${ADAPTER}" -a ! -z "${READDIR}" ]; then
	if [ -f "$ADAPTER" ]; then
        ## USE FASTQ-MCF ##
        fastq-mcf -o $READDIR/clipped_$ID.fastq -l $MINLEN -q $MINQUAL -w 4 -x 10 -t 0 $ADAPTER $READ1

        # trim first few nucleotides
        TRIMLEN=$(($TRIMLEN + 1))
        fastx_trimmer -Q 33 -f $TRIMLEN -i $READDIR/clipped_$ID.fastq -m $MINLEN -o $READDIR/clipped_$ID.fastq.tmp
        mv $READDIR/clipped_$ID.fastq.tmp $READDIR/clipped_$ID.fastq

        # rename trimmed file from tmp to original
		if [ ! -z "$FASTA" ]; then
            fastq_to_fasta -Q 33 -n -i $READDIR/clipped_$ID.fastq -o $READDIR/clipped_$ID.fasta 
        fi

		## USE TRIMMOMATIC ##
		#eval $trimSE -phred33 \"$READ1\" $READDIR/$ID"".fastq.clipped.gz ILLUMINACLIP:$ADAPTER:2:30:10 LEADING:20 TRAILING:20 MINLEN:18
		## convert fastq to fasta
		#zless "$READDIR/$ID"".fastq.clipped.gz" | fastq_to_fasta -n | perl -ane 'chomp($_); if($_=~/^\>/) { if(defined($header)) { print "$header#".length($seq)."\n$seq\n"; } $header=$_; $seq=(); } else { $seq.=$_; } END { print "$header|".length($seq)."\n$seq\n"; }' > "$READDIR/$ID"".fasta"

		## USE FASTX ##
		#COMMAND_CLIP=`zgrep -v "^>" $ADAPTER | perl -ane 'chomp($_); $_=~s/[^a-zA-Z]+//g; $command.=sprintf("fastx_clipper -a %s -l 18 -Q 33 -n |", $_); END { $command=~s/\|$//g; print $command."\n"; }'`;
		#if [ ! -z "$FASTA" ]; then
		#	zless $READ1 | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.fastq.clipped
		#	zless $READDIR/$ID.fastq.clipped | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#else
		#	zless $READ1 | $COMMAND_CLIP | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		#fi
	else
		#zless $READ1 | fastx_clipper -a $ADAPTER -l 18 -Q 33 | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
		## use fastx_clipper with '-n' parameter, if fasta file is empty
		if [ ! -s "$READDIR/$ID.fasta" ]; then
			if [ ! -z "$FASTA" ]; then
				## remove beginning of the read until the 'N'
				zless $READ1 | fastx_clipper -a $ADAPTER -l 18 -Q 33 -n | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 > $READDIR/$ID.fastq.clipped
				zless $READDIR/$ID.fastq.clipped | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta
			else
				## remove beginning of the read until the 'N'
				zless $READ1 | fastx_clipper -a $ADAPTER -l 18 -Q 33 -n | fastq_quality_trimmer -t 20 -l 18 -Q 33 | fastx_artifacts_filter -Q 33 | fastq_to_fasta -Q 33 -n | perl -ane 'if($_=~/^\>/) { $header=$_; } else { $_=~s/^.+N//g; if(length($_)>=18) { print "$header$_"; } }' | fastx_collapser -Q 33 | perl -ane 'if($_=~/^>/) { $_=~s/\>//g; @t=split(/\-/,$_); print ">'$ID'_$t[0]|$t[1]"; } else { print $_; }' > $READDIR/$ID.fasta

			fi
		fi
	fi
fi

## map reads
if [ ! -z "${MAPDIR}" ]; then
	#echo "$FASTAFILE $GENOMEINDEX $READDIR $ID"; exit;
	if [ -e "$FASTAFILE" -a -e "$GENOMEINDEX" -a -e "$READDIR/$ID.fasta" -o -e "$READDIR/clipped_$ID.fastq" ]; then
		#tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $READDIR/clipped_$ID.fastq

        ## compute mapping statistics
        #samtools index $MAPDIR/$ID/accepted_hits.bam $MAPDIR/$ID/accepted_hits.bai && samtools idxstats $MAPDIR/$ID/accepted_hits.bam > $MAPDIR/$ID/accepted_MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID/accepted_MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt &
        #samtools index $MAPDIR/$ID/unmapped.bam $MAPDIR/$ID/unmapped.bai && samtools idxstats $MAPDIR/$ID/unmapped.bam > $MAPDIR/$ID/unmapped_MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID/unmapped_MappingStatistics.txt >> $MAPDIR/concatenated_unmapped_MappingStatistics.txt &

        ## create bigwig files for viualization at the UCSC genome browser
        bedtools bamtobed -i $MAPDIR/$ID/accepted_hits.bam -bed12 | grep '^[1-9XY]' | awk '{print "chr"$0}' > $MAPDIR/$ID/accepted_hits_corrected.bed && bedtools genomecov -bg -i $MAPDIR/$ID/accepted_hits_corrected.bed -g $CHRSIZE -split > $MAPDIR/$ID/accepted_hits.bedGraph && bedGraphToBigWig $MAPDIR/$ID/accepted_hits.bedGraph $CHRSIZE $UCSCDIR/$ID.bw && rm $MAPDIR/$ID/accepted_hits.bedGraph

		#mkdir $MAPDIR/$ID
		#bowtie2 -p 8 -x $GENOMEINDEX -U $READDIR/$ID".fasta" -f -S $MAPDIR/$ID/$ID.sam 

        ## MEDIP-seq
		#segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta" -V -A 95 > $MAPDIR/$ID.sam
		#segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta1" -V -A 95 > $MAPDIR/$ID.sam1

        ## small RNA-seq
		#segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta" -V -A 85 > $MAPDIR/$ID.sam
		#bam2bed.pl -i $MAPDIR/$ID.sam -s -o $MAPDIR/$ID.bed

	fi
fi

echo "Done"
