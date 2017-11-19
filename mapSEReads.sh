#PBS -l nodes=1:ppn=4

## DEPENDENCIES
MAPDIR="."
PROCESSORS=1
GENOME="mm9"
ALNMODE="--sensitive"
TRIM5=0
TRIM3=0

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
    echo " -p <int>    [number of processors (default: 1)]"
    echo " -d <string> [identifier for output BAM file (default: same as fastq file)]"
    echo " -s          [perform alignment accommodating for splice junctions using tophat2]"
    echo "             [default is to use bowtie2]"
    echo " -r          [map reads for repeat analysis using bowtie (output multimapped reads)]"
    echo " -R          [map reads for repeat analysis using bowtie2]"
    echo "[OPTIONS: bowtie2]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -e          [extend 3' end of reads in output bigWig files]"
    echo " -k <int>    [instead of reporting best alignment, report input number of alignments per read]"
    echo " -q <string> [end-to-end: --very-fast, --fast, --sensitive, --very-sensitive (default: --sensitive)]"
    echo "             [local: --very-fast-local, --fast-local, --sensitive-local, --very-sensitive-local (default: --sesitive-local)"
    echo " -l          [local alignment; ends might be soft clipped]"
    echo " -f <int>    [trim <int> bases from 5'/left end of reads (default: 0)]"
    echo " -t <int>    [trim <int> bases from 3'/right end of reads (default: 0)]"
    echo " -L <int>    [length of seed substring; must be >3, <32 (default: 22)]"
    echo " -I <string> [interval between seed substrings w/r/t read len (default: S,1,1.15)]"
    echo " -D <int>    [give up extending after <int> failed extends in a row (default: 15)]"
    echo " -R <int>    [for reads w/ repetitive seeds, try <int> sets of seeds (default: 2)]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:m:g:p:d:srRucek:q:lf:t:L:I:D:R:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		m) MAPDIR=$OPTARG;;
		g) GENOME=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
        d) ID=$OPTARG;;
        s) SPLICE=1;;
        r) REPENRICH=1;;
        R) REPEATS=1;;
        u) UNIQUE=1;;
        c) SCALE=1;;
        e) EXTEND=1;;
        k) ALNCOUNT=$OPTARG;;
        q) ALNMODE=$OPTARG;;
        l) LOCAL=1;;
        f) TRIM5=$OPTARG;;
        t) TRIM3=$OPTARG;;
        L) SEED=$OPTARG;;
        I) INTERVAL=$OPTARG;;
        D) GIVEUP=$OPTARG;;
        R) TRIES=$OPTARG;;
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
elif [ "$GENOME" == "hg19_ifn" ]; then
    GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/interferon_genes/interferon"
    FASTAFILE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/interferon_genes/interferon"
    CHRSIZE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/ChromInfoRef.txt"
else
    echo "Presently the program only support analysis for mm9 or hg19"
    echo
    usage
fi
echo done

## retrieve file name
if [ -z "$ID" ]; then
    ID=`echo $FASTQ | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
fi
FASTQ=$(echo $FASTQ | sed 's/\,/ /g')
READLENGTH=`zless $FASTQ | head -n 2 | tail -n 1 | perl -ane '$len=length($_)-1; print $len;'`;
#echo -e "$ID\t$READLENGTH"; exit;

## read arguments
ARGS=""
if [ ! -z "$ALNCOUNT" ]; then
    ARGS="$ARGS -k $ALNCOUNT";
fi

if [ ! -z "$LOCAL" ]; then
    ARGS="$ARGS --local";
fi

if [ ! -z "$SEED" ]; then
    ARGS="$ARGS -L $SEED";
fi

if [ ! -z "$INTERVAL" ]; then
    ARGS="$ARGS -i $INTERVAL";
fi

if [ ! -z "$GIVEUP" ]; then
    ARGS="$ARGS -D $GIVEUP";
fi

if [ ! -z "$TRIES" ]; then
    ARGS="$ARGS -R $TRIES";
fi

## map reads
echo "Map for $ID... " >$MAPDIR/$ID.mapStat
#echo "$FASTAFILE $GENOMEINDEX $READDIR $ID"; exit;

## start analysis
if [ ! -z "$SPLICE" ]; then
    echo "Command used: tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ" >>$MAPDIR/$ID.mapStat
    tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ

    ## compute mapping statistics
    #samtools index $MAPDIR/$ID"_accepted_hits.bam" $MAPDIR/$ID"_accepted_hits.bai" && samtools idxstats $MAPDIR/$ID"_accepted_hits.bam" > $MAPDIR/$ID"_accepted_MappingStatistics.txt" && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID"_accepted_MappingStatistics.txt" >> $MAPDIR/concatenated_accepted_MappingStatistics.txt &
    #samtools index $MAPDIR/$ID"_unmapped.bam" $MAPDIR/$ID"_unmapped.bai" && samtools idxstats $MAPDIR/$ID"_unmapped.bam" > $MAPDIR/$ID"_unmapped_MappingStatistics.txt" && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID"_unmapped_MappingStatistics.txt" >> $MAPDIR/concatenated_unmapped_MappingStatistics.txt &
    mv $MAPDIR/$ID/accepted_hits.bam $MAPDIR/$ID.bam
    samtools index $MAPDIR/$ID.bam
    zless $MAPDIR/$ID/align_summary.txt >> $MAPDIR/$ID.mapStat

    ## create bigwig files for viualization at the UCSC genome browser
    bedtools bamtobed -i $MAPDIR/$ID.bam -bed12 | grep '^[1-9XY]' | awk '{print "chr"$0}' > $MAPDIR/$ID/accepted_hits_corrected.bed && bedtools genomecov -bg -i $MAPDIR/$ID/accepted_hits_corrected.bed -g $CHRSIZE -split > $MAPDIR/$ID/accepted_hits.bedGraph && bedGraphToBigWig $MAPDIR/$ID/accepted_hits.bedGraph $CHRSIZE $MAPDIR/$ID.bw && rm $MAPDIR/$ID/accepted_hits.bedGraph
elif [ ! -z "$REPENRICH" ]; then
    if [ ! -d "$MAPDIR" ]; then
        mkdir $MAPDIR/
    fi

    bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID"_multimap.fastq" $FASTQ $MAPDIR/$ID"_unique.sam" 2>>$MAPDIR/$ID.mapStat
    samtools view -bS $MAPDIR/$ID"_unique.sam" > $MAPDIR/$ID"_unique.bam"
    samtools sort $MAPDIR/$ID"_unique.bam" -o $MAPDIR/$ID"_unique_sorted.bam"
    mv $MAPDIR/$ID"_unique_sorted.bam" $MAPDIR/$ID"_unique.bam"
    samtools index $MAPDIR/$ID"_unique.bam"
    #samtools idxstats $MAPDIR/$ID"_unique.bam" > $MAPDIR/$ID.MappingStatistics.txt
    rm $MAPDIR/$ID"_unique.sam"
elif [ ! -z "$REPEATS" ]; then
    if [ ! -d "$MAPDIR" ]; then
        mkdir $MAPDIR/
    fi

    ## inspired from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1375157
    ## command check
    echo "Command used: zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000" >>$MAPDIR/$ID.mapStat

    zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort - -o $MAPDIR/$ID.bam 

    samtools index $MAPDIR/$ID.bam

    ## create bigwig files for visualization at the UCSC genome browser
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
else
    if [ ! -d "$MAPDIR" ]; then
        mkdir $MAPDIR/
    fi

<<"COMMENT"
COMMENT
    ## command check
    echo "Command used: zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS" >>$MAPDIR/$ID.mapStat

    if [ ! -z "$UNIQUE" ]; then
            zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | grep -v XS: | samtools view -S -b - | samtools sort - -o $MAPDIR/$ID.bam
    else
        zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort - -o $MAPDIR/$ID.bam
    fi

    ## compute mapping statistics
    ## idxstat format: The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
    samtools index $MAPDIR/$ID.bam
    #samtools index $MAPDIR/$ID.bam && samtools idxstats $MAPDIR/$ID.bam > $MAPDIR/$ID.MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID.MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt

    ## create bigwig files for visualization at the UCSC genome browser
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
