#PBS -l nodes=1:ppn=4

## DEPENDENCIES
MAPDIR="."
PROCESSORS=1
GENOME="mm9"
ALNMODE="--sensitive"
TRIM5=0
TRIM3=0
KALLISTO_FL=200
KALLISTO_SD=30

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
	echo " -m <dir>    [output directory to store mapped reads (default: .)]"
	echo " -g <string> [genome (default: mm9)]"
    echo "             [mm9, mm10, hg19, hg38, dm6, ERCC]"
    echo " -p <int>    [number of processors (default: 1)]"
    echo " -d <string> [identifier for output BAM file (default: same as fastq file)]"
    echo " -a <dir>    [directory to keep adapter trimmed fastq files]"
    echo "             [required to auto-trim adapter sequences]"
    echo " -z <dir>    [directory to keep spike-in normalization files]"
    echo "             [required to auto-determine scaling factor using splike-in data]"
    echo " -Z          [only bam to bw coversion; bam file exists]"
    echo "[OPTIONS: repeats]"
    echo " -r          [map reads for repeat analysis using bowtie (output multimapped reads)]"
    echo " -R          [map reads for repeat analysis using bowtie2]"
    echo "[OPTIONS: tophat2]"
    echo " -s          [perform alignment accommodating for splice junctions using tophat2]"
    echo "[OPTIONS: STAR]"
    echo " -S          [perform alignment accommodating for splice junctions using STAR]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -f <int>    [trim <int> bases from 5'/left end of reads (default: 0)]"
    echo " -t <int>    [trim <int> bases from 3'/right end of reads (default: 0)]"
    echo "[OPTIONS: bowtie2 (default)]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files (overridden if -z)]"
    echo " -C          [scale the read coverage to 1x in output bigWig files (overridden if -z)]"
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
    echo " -E <int>    [for reads w/ repetitive seeds, try <int> sets of seeds (default: 2)]"
    echo "[OPTIONS: Kallisto]"
    echo " -K          [perform alignment using kallisto]"
    echo " -T <int>    [average fragment length (default: 200)]"
    echo " -N <int>    [standard deviation of fragment length (default: 30)]"
    echo "[NOTE: splike-in scaling formula]"
    echo "             [https://www.sciencedirect.com/science/article/pii/S2211124714008729#app3]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:m:g:p:d:a:z:ZrRsSucCek:q:lf:t:L:I:D:E:KT:N:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		m) MAPDIR=$OPTARG;;
		g) GENOME=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
        d) ID=$OPTARG;;
        a) TRIM_FASTQ_DIR=$OPTARG;;
        z) SPIKEIN=$OPTARG;;
        Z) BAMTOBW=1;;
        r) REPENRICH=1;;
        R) REPEATS=1;;
        s) TOPHAT2=1;;
        S) STAR=1;;
        u) UNIQUE=1;;
        c) SCALE=1;;
        C) SCALE1x=1;;
        e) EXTEND=1;;
        k) ALNCOUNT=$OPTARG;;
        q) ALNMODE=$OPTARG;;
        l) LOCAL=1;;
        f) TRIM5=$OPTARG;;
        t) TRIM3=$OPTARG;;
        L) SEED=$OPTARG;;
        I) INTERVAL=$OPTARG;;
        D) GIVEUP=$OPTARG;;
        E) TRIES=$OPTARG;;
        K) KALLISTO=1;;
        T) KALLISTO_FL=$OPTARG;;
        N) KALLISTO_SD=$OPTARG;;
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
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/kallisto/Mus_musculus.NCBIM37.67.cdna.all.idx"
    else
        ## bowtie2 (*bt2)
        #GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie2/Bowtie2IndexWithAbundance"
        ## bowtie2 (*bt2 - with chromosome - no scaffolds)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie2_chr_noscaffold/Bowtie2IndexWithAbundance"
        FASTAFILE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/TopHatTranscriptomeIndex_with_chr/genes_without_mt"
        CHRSIZE="/home/pundhir/software/RNAPipe/data/Mus_musculus/Ensembl/NCBIM37/ChromInfoRef.txt"
    fi
elif [ "$GENOME" == "mm10" ]; then
    if [ ! -z "$REPENRICH" ]; then
        ## tophat (bowtie1 - *ebwt)
        GENOMEINDEX=""
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/mm10/kallisto/Mus_musculus.GRCm38.cdna.all.idx"
    else
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Mus_musculus/mm10/bowtie2/Bowtie2IndexWithAbundance"
        FASTAFILE=""
        CHRSIZE="/home/pundhir/project/genome_annotations/mouse.mm10.genome"
    fi
elif [ "$GENOME" == "hg19" ]; then
    if [ ! -z "$REPENRICH" ]; then
        ## tophat (bowtie1 - *ebwt)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie/genome_and_Abundant"
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37"
    else
        ## bowtie2 (*bt2)
        #GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie2/genome_and_Abundant"
        ## bowtie2 (*bt2 - with chromosome - no scaffolds)
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie2_chr_noscaffold/genome_and_Abundant"
        FASTAFILE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/TopHatTranscriptomeIndex_with_chr/genes_without_mt"
        CHRSIZE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/ChromInfoRef.txt"
    fi
elif [ "$GENOME" == "hg38" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/hg38/kallisto/Homo_sapiens.GRCh38.cdna.all.idx"
    else
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/hg38/bowtie2/Bowtie2IndexWithAbundance"
        FASTAFILE=""
        CHRSIZE="/home/pundhir/project/genome_annotations/human.hg38.genome"
    fi
elif [ "$GENOME" == "hg19_ifn" ]; then
    GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Homo_sapiens/interferon_genes/interferon"
    FASTAFILE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/interferon_genes/interferon"
    CHRSIZE="/home/pundhir/software/RNAPipe/data/Homo_sapiens/Ensembl/GRCh37/ChromInfoRef.txt"
elif [ "$GENOME" == "dm6" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Dro_melanogaster/dm6/kallisto/"
    else
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/Dro_melanogaster/dm6/bowtie2/Bowtie2IndexWithAbundance"
        FASTAFILE=""
        CHRSIZE="/home/pundhir/project/genome_annotations/drosophila.dm6.genome"
    fi
elif [ "$GENOME" == "ERCC" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/ERCC/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/ERCC/kallisto/"
    else
        GENOMEINDEX="/home/pundhir/software/RNAPipe/data/ERCC/bowtie2/Bowtie2IndexWithAbundance"
        FASTAFILE=""
        CHRSIZE=""
    fi
else
    echo "Presently the program only support analysis for mm9, mm10, hg19, hg38, dm6 or ERCC"
    echo
    usage
fi
echo done

## auto-trim adapter sequences, if required
if [ ! -z "$TRIM_FASTQ_DIR" -a -z "$BAMTOBW" ]; then
    ## parse input fastq files in an array
    oIFS=$IFS
    IFS=","
    FASTQFILES=($FASTQ)
    FASTQFILES_COUNT=${#FASTQFILES[@]}
    IFS=$oIFS

    FASTQ=""
    for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
        TEMPID=`echo ${FASTQFILES[$i]} | perl -an -F'/\,/' -e '$TEMPID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $TEMPID.=$_."_"; } $TEMPID=~s/\_$//g; print "$TEMPID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
        if [ ! -f "$TRIM_FASTQ_DIR/$TEMPID.clipped.fastq.gz" ]; then
            trimAdapters.sh -i ${FASTQFILES[$i]} -r $TRIM_FASTQ_DIR -a auto
        fi
        FASTQ="$FASTQ,$TRIM_FASTQ_DIR/$TEMPID.clipped.fastq.gz"
    done
    FASTQ=$(echo $FASTQ | perl -ane '$_=~s/^\,//g; print $_;')
fi

## retrieve file name
if [ -z "$ID" ]; then
    ID=`echo $FASTQ | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
fi
FASTQ=$(echo $FASTQ | sed 's/\,/ /g')
READLENGTH=`zless $FASTQ | head -n 2 | tail -n 1 | perl -ane '$len=length($_)-1; print $len;'`;
#echo -e "$ID\t$READLENGTH"; exit;

## splike-in normalization
if [ ! -z "$SPIKEIN" ]; then
    ARGS=""
    if [ ! -z "$UNIQUE" ]; then
        ARGS="$ARGS -u";
    fi
    if [ ! -z "$EXTEND" ]; then
        ARGS="$ARGS -e";
    fi

    if [ -z "$SPIKEIN/$ID.mapStat" ]; then
        mapSEReads.sh -i $FASTQ -g dm6 -m $SPIKEIN -p $PROCESSORS -d $ID $ARGS
    fi

    if [ -z "$UNIQUE" ]; then
        SCALE_SPIKEIN=$(zless $SPIKEIN/$ID.mapStat | grep 'aligned exactly\|aligned >1' | perl -ane '$sum+=$F[0]; END { print "$sum\n"; }')
    else
        SCALE_SPIKEIN=$(zless $SPIKEIN/$ID.mapStat | grep 'aligned exactly' | perl -ane '$sum+=$F[0]; END { print "$sum\n"; }')
    fi

    SCALE_SPIKEIN=$(echo $SCALE_SPIKEIN | perl -ane 'printf("%0.6f", 1000000/$_);');
fi
#echo $SCALE_SPIKEIN; exit

## read arguments
if [ ! -z "$STAR" ]; then
    ARGS=""
    if [ ! -z "$SCALE" ]; then
        ARGS="$ARGS --outWigNorm RPM";
    else
        ARGS="$ARGS --outWigNorm None";
    fi
else
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
fi

#echo "$FASTAFILE $GENOMEINDEX $READDIR $ID"; exit;

## start analysis
if [ ! -z "$STAR" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        FASTQ=$(echo $FASTQ | sed 's/ /\,/g')

        echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS" >> $MAPDIR/$ID.mapStat

        STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS
        mv $MAPDIR/$ID"Aligned.sortedByCoord.out.bam" $MAPDIR/$ID.bam
        zless $MAPDIR/$ID"Log.final.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.progress.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.out" > $MAPDIR/$ID.log
        zless $MAPDIR/$ID"SJ.out.tab" > $MAPDIR/$ID.SJ
        samtools index $MAPDIR/$ID.bam
        rm $MAPDIR/$ID"Log.final.out" $MAPDIR/$ID"Log.progress.out" $MAPDIR/$ID"Log.out" $MAPDIR/$ID"SJ.out.tab"

        if [ ! -z "$UNIQUE" ]; then
            mv $MAPDIR/$ID"Signal.Unique.str1.out.bg" $MAPDIR/$ID.bw
            rm $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg" 
        else
            mv $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg" $MAPDIR/$ID.bw
            rm $MAPDIR/$ID"Signal.Unique.str1.out.bg"
        fi
    fi
elif [ ! -z "$TOPHAT2" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        echo "Command used: tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ" >>$MAPDIR/$ID.mapStat
        tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ

        ## compute mapping statistics
        mv $MAPDIR/$ID/accepted_hits.bam $MAPDIR/$ID.bam
        samtools index $MAPDIR/$ID.bam
        zless $MAPDIR/$ID/align_summary.txt >> $MAPDIR/$ID.mapStat
    fi

    ## create bigwig files for viualization at the UCSC genome browser
    bedtools bamtobed -i $MAPDIR/$ID.bam -bed12 | grep '^[1-9XY]' | awk '{print "chr"$0}' > $MAPDIR/$ID/accepted_hits_corrected.bed && bedtools genomecov -bg -i $MAPDIR/$ID/accepted_hits_corrected.bed -g $CHRSIZE -split > $MAPDIR/$ID/accepted_hits.bedGraph && bedGraphToBigWig $MAPDIR/$ID/accepted_hits.bedGraph $CHRSIZE $MAPDIR/$ID.bw && rm $MAPDIR/$ID/accepted_hits.bedGraph
elif [ ! -z "$REPENRICH" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        echo "Command used: bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID'_multimap.fastq' $FASTQ $MAPDIR/$ID'_unique.sam'" >> $MAPDIR/$ID.mapStat

        bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID"_multimap.fastq" $FASTQ $MAPDIR/$ID"_unique.sam" 2>>$MAPDIR/$ID.mapStat
        samtools view -bS $MAPDIR/$ID"_unique.sam" > $MAPDIR/$ID"_unique.bam"
        samtools sort $MAPDIR/$ID"_unique.bam" -o $MAPDIR/$ID"_unique_sorted.bam"
        mv $MAPDIR/$ID"_unique_sorted.bam" $MAPDIR/$ID"_unique.bam"
        samtools index $MAPDIR/$ID"_unique.bam"
        #rm $MAPDIR/$ID"_unique.sam"
    fi
elif [ ! -z "$REPEATS" ]; then
    if [ -z "$BAMTOBW" ]; then
        ## inspired from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1375157
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        ## command check
        echo "Command used: zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000" >>$MAPDIR/$ID.mapStat

        zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort -m 1500M - -o $MAPDIR/$ID.bam 

        samtools index $MAPDIR/$ID.bam
    fi
elif [ ! -z "$KALLISTO" ]; then
    if [ -z "$BAMTOBW" ]; then
        ## customize output directory in case it is not provided
        if [ "$MAPDIR" -eq "." ]; then
            MAPDIR=$ID
            mkdir $MAPDIR
        fi
        kallisto quant -i $GENOMEINDEX -o $MAPDIR -b 100 --single --bias -l $KALLISTO_FL -s $KALLISTO_SD -t $PROCESSORS $FASTQ
    fi
else
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        ## command check
        echo "Command used: zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS" >>$MAPDIR/$ID.mapStat

        if [ ! -z "$UNIQUE" ]; then
                zless $FASTQ  | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | grep -v XS: | samtools view -S -b - | samtools sort -m 1500M - -o $MAPDIR/$ID.bam
        else
            zless $FASTQ | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort -m 1500M - -o $MAPDIR/$ID.bam
        fi

        ## compute mapping statistics
        ## idxstat format: The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
        #samtools index $MAPDIR/$ID.bam && samtools idxstats $MAPDIR/$ID.bam > $MAPDIR/$ID.MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID.MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt
        samtools index $MAPDIR/$ID.bam
    fi
fi 

## create bigwig files for visualization at the UCSC genome browser
if [ ! -z "$REPENRICH" -o ! -z "$REPEATS" -o -z "$STAR" -a -z "$TOPHAT2" -a -z "$KALLISTO" ]; then
    if [ ! -z "$SPIKEIN" ]; then
        echo -e "Using spike-in scale, $SCALE_SPIKEIN to normalize bigWig files.. "
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -z $SCALE_SPIKEIN -p $PROCESSORS
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -z $SCALE_SPIKEIN -p $PROCESSORS
        fi
        echo "done"
    elif [ ! -z "$SCALE" ]; then
        echo -e "Using RPKM scaling to normalize bigWig files.. "
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -c -p $PROCESSORS
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -c -p $PROCESSORS
        fi
        echo "done"
    elif [ ! -z "$SCALE1x" ]; then
        echo -e "Using 1x scaling to normalize bigWig files.. "
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -C -p $PROCESSORS
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -C -p $PROCESSORS
        fi
        echo "done"
    else
        echo -e "Using no scaling to normalize bigWig files.. "
        if [ ! -z "$EXTEND" ]; then
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -p $PROCESSORS
        else
            bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -p $PROCESSORS
        fi
        echo "done"
    fi
fi
echo "done"
