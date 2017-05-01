#PBS -l nodes=1:ppn=4

## DEPENDENCIES
MAPDIR="."
PROCESSORS=1
GENOME="mm9"

#### usage ####
usage() {
	echo Program: "mapSEReads.sh (map single end reads)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@bric.ku.dk
	echo "Usage: mapSEReads.sh -i <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input fastq file(s) with single end reads]"
    echo "             [if multiple separate them by a comma]"
	echo "[OPTIONS]"
	echo " -m <dir>    [directory to store mapped reads (default: .)]"
	echo " -g <string> [genome (default: mm9)]"
    echo "             [mm9 or hg19]"
    echo " -s          [perform alignment accommodating for splice junctions using tophat2]"
    echo "             [default is to use bowtie2]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -e          [extend 3' end of reads in output bigWig files]"
    echo " -p <int>    [number of processors (default: 1)]"
    echo " -r          [map reads for repeat analysis using RepEnrich]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:m:g:sucep:rh ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		m) MAPDIR=$OPTARG;;
		g) GENOME=$OPTARG;;
        s) SPLICE=1;;
        u) UNIQUE=1;;
        c) SCALE=1;;
        e) EXTEND=1;;
        p) PROCESSORS=$OPTARG;;
        r) REPENRICH=1;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$FASTQ" -o "$HELP" ]; then
	usage
fi

## create appropriate directory structure
echo -n "Create appropriate directory structure... "
if [ ! -d "$MAPDIR" ]; then
    mkdir $MAPDIR
fi
echo done

echo -n "Populating files based on input genome, $GENOME (`date`).. "
if [ "$GENOME" == "mm9" ]; then
    if [ ! -z "$REPENRICH" ]; then
        ## tophat (bowtie1 - *ebwt)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie/Bowtie2IndexWithAbundance"
    else
        ## bowtie2 (*bt2)
        #GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie2/Bowtie2IndexWithAbundance"
        ## bowtie2 (*bt2 - with chromosome)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie2_chr_noscaffold/Bowtie2IndexWithAbundance"
        FASTAFILE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/TopHatTranscriptomeIndex/bowtie2/genes_without_mt"
        CHRSIZE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/ChromInfoRef.txt"
    fi
elif [ "$GENOME" == "hg19" ]; then
    if [ ! -z "$REPENRICH" ]; then
        ## tophat (bowtie1 - *ebwt)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie/genome_and_Abundant"
    else
        ## bowtie2 (*bt2)
        #GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie2/genome_and_Abundant"
        ## bowtie2 (*bt2 - with chromosome)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie2_chr/genome_and_Abundant"
        FASTAFILE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/TopHatTranscriptomeIndex/bowtie2/genes_without_mt"
        CHRSIZE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/ChromInfoRef.txt"
    fi
else
    echo "Presently the program only support analysis for mm9 or hg19"
    echo
    usage
fi
echo done

## retrieve file name
ID=`echo $FASTQ | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
FASTQ=$(echo $FASTQ | sed 's/\,/ /g')
READLENGTH=`zless $FASTQ | head -n 2 | tail -n 1 | perl -ane '$len=length($_)-1; print $len;'`;
#echo -e "$ID\t$READLENGTH"; exit;

## map reads
echo "Map for $ID... "
#echo "$FASTAFILE $GENOMEINDEX $READDIR $ID"; exit;
if [ ! -z "$SPLICE" ]; then
    tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ

    ## compute mapping statistics
    samtools index $MAPDIR/$ID/accepted_hits.bam $MAPDIR/$ID/accepted_hits.bai && samtools idxstats $MAPDIR/$ID/accepted_hits.bam > $MAPDIR/$ID/accepted_MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID/accepted_MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt &
    samtools index $MAPDIR/$ID/unmapped.bam $MAPDIR/$ID/unmapped.bai && samtools idxstats $MAPDIR/$ID/unmapped.bam > $MAPDIR/$ID/unmapped_MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID/unmapped_MappingStatistics.txt >> $MAPDIR/concatenated_unmapped_MappingStatistics.txt &

    ## create bigwig files for viualization at the UCSC genome browser
    bedtools bamtobed -i $MAPDIR/$ID/accepted_hits.bam -bed12 | grep '^[1-9XY]' | awk '{print "chr"$0}' > $MAPDIR/$ID/accepted_hits_corrected.bed && bedtools genomecov -bg -i $MAPDIR/$ID/accepted_hits_corrected.bed -g $CHRSIZE -split > $MAPDIR/$ID/accepted_hits.bedGraph && bedGraphToBigWig $MAPDIR/$ID/accepted_hits.bedGraph $CHRSIZE $MAPDIR/$ID.bw && rm $MAPDIR/$ID/accepted_hits.bedGraph
elif [ ! -z "$REPENRICH" ]; then
    if [ ! -d "$MAPDIR" ]; then
        mkdir $MAPDIR/
    fi

    bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID"_multimap.fastq" $FASTQ $MAPDIR/$ID"_unique.sam"
    samtools view -bS $MAPDIR/$ID"_unique.sam" > $MAPDIR/$ID"_unique.bam"
    samtools sort $MAPDIR/$ID"_unique.bam" -o $MAPDIR/$ID"_unique_sorted.bam"
    mv $MAPDIR/$ID"_unique_sorted.bam" $MAPDIR/$ID"_unique.bam"
    samtools index $MAPDIR/$ID"_unique.bam"
    rm $MAPDIR/$ID"_unique.sam"
else
    if [ ! -d "$MAPDIR" ]; then
        mkdir $MAPDIR/
    fi

<<"COMMENT"
COMMENT
    if [ ! -z "$UNIQUE" ]; then
        zcat -f $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - | grep -v XS: | samtools view -S -b - | samtools sort - -o $MAPDIR/$ID.bam &>$MAPDIR/$ID.log
    else
        zcat -f $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - | samtools view -S -b - | samtools sort - -o $MAPDIR/$ID.bam &>$MAPDIR/$ID.log
    fi

    ## compute mapping statistics
    ## idxstat format: The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
    samtools index $MAPDIR/$ID.bam && samtools idxstats $MAPDIR/$ID.bam > $MAPDIR/$ID.MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID.MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt &

    ## create bigwig files for viualization at the UCSC genome browser
    if [ ! -z "$SCALE" ]; then
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -s
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -s
        fi
    else
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME
        fi
    fi

    ## MEDIP-seq
    #segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta" -V -A 95 > $MAPDIR/$ID.sam
    #segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta1" -V -A 95 > $MAPDIR/$ID.sam1

    ## small RNA-seq
    #segemehl.x -s --minsize 18 -t 8 -d $FASTAFILE -i $GENOMEINDEX -q $READDIR/$ID".fasta" -V -A 85 > $MAPDIR/$ID.sam
    #bam2bed.pl -i $MAPDIR/$ID.sam -s -o $MAPDIR/$ID.bed
fi 

echo "done"
