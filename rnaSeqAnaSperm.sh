#!/bin/bash
#PBS -l nodes=1:ppn=4

ORG="hsa"
ORGDIR="sperm"

MAP=/home/projects/deepalign/data/maps/$ORG/$ORGDIR
CLUSTER=/home/projects/deepalign/data/clusters/$ORG/$ORGDIR
UCSC=/home/projects/deepalign/data/ucsc/$ORG/$ORGDIR
BLOCKSTAT=/home/projects/deepalign/data/compBlockStat/$ORG/$ORGDIR
ANNOTATION=/home/projects/deepalign/data/annotations/$ORG
SOFTWARE=/home/users/sachin/software/myScripts

FASTQ=("$MAP/spermInfertile1.bed.gz" , "$MAP/spermInfertile2.bed.gz")

if [ ! -d "$MAP" -o ! -d "$CLUSTER" -o ! -d "$UCSC" -o ! -d "$BLOCKSTAT" -o ! -d "$ANNOTATION" -o ! -d "$SOFTWARE" ];
then
	echo "Cannot locate one or more input directories"
	exit
fi

if [ -f $CLUSTER/bgAnnoDist.out ];
then
	rm $CLUSTER/bgAnnoDist.out
fi

#for fastq in $MAP/*.bed.gz
for fastq in "${FASTQ[@]}"
do
	cd /home/projects/deepalign/data/reads/hsa/humanBodyMap/tags
	id=`echo $fastq | perl -ane '@L=split(/\//,$F[0]);foreach(@L){if(/bed\.gz/){($id)=split(/\./,$_);print "$id\n";}}'`

	#### Step-1: sort bed file by (strand, chr, start and end)
	if [ ! -f $MAP/$id.bed.gz ];
	then
		echo "Cannot locate $MAP/$id.bed.gz"
		exit
	fi
	#zcat $MAP/$id.bed.gz | perl -an -F'/\s+/' -e 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t1\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[10]\t$F[11]\n";' | sort -k 6,6 -k 1,1 -k 2n,2 -k 3n,3 -o $MAP/$id.bed.sorted

	#### Step-2: define block and block groups for each chromosome using blockbuster ####
	if [ -f $MAP/$id.bed.sorted.gz ];
	then
		gunzip $MAP/$id.bed.sorted.gz
	fi
	if [ ! -f $MAP/$id.bed.sorted ];
	then
		echo "Cannot locate $MAP/$id.bed.sorted"
		exit
	fi
	
	#CHR=(`cut -f 1 $MAP/$id.bed.sorted | grep -v "_" | uniq | sort | uniq`)

	## remove output file, if exists
	if [ -f $CLUSTER/$id.map.clusters.gz ];
	then
		rm $CLUSTER/$id.map.clusters.gz
	fi
	if [ -f $CLUSTER/$id.map.clusters ];
	then
		rm $CLUSTER/$id.map.clusters
	fi

	#for chr in "${CHR[@]}"
	#do
		### humanBodyMap
		#blockbuster.x -format 2 -minBlockHeight 10 -minClusterHeight 50 -scale 0.5 -chr $chr -print 2 $MAP/$id.map.sorted >> $CLUSTER/$id.map.clusters
		### emtab305 and egeod28123
		#blockbuster.x -format 2 -minBlockHeight 2 -minClusterHeight 10 -scale 0.5 -chr $chr -print 2 $MAP/$id.map.sorted >> $CLUSTER/$id.map.clusters
		#blockbuster.x -format 1 -minBlockHeight 10 -minClusterHeight 50 -scale 0.5 -chr $chr -print 2 $MAP/$id.bed.sorted >> $CLUSTER/$id.map.clusters
	#done
	
	cat $MAP/$id.bed.sorted | perl -an -F'/\s+/' -e 'BEGIN{$end=100000; open(OUTFILE, ">'$MAP/$id.tmp'"); } if(($F[1] - $end) < 200) { print OUTFILE $_; $end=$F[2]; } else { close OUTFILE; system("blockbuster.x -format 1 -minBlockHeight 10 -minClusterHeight 50 -scale 0.5 -print 2 '$MAP/$id.tmp' >>'$CLUSTER/$id.map.clusters'"); open(OUTFILE, ">'$MAP/$id.tmp'"); $end=$F[2]; print OUTFILE $_; }'

	## validate block groups defined by blockbuster
	validateBlockbuster.pl -c $CLUSTER/$id.map.clusters -d 60 > $CLUSTER/$id.map.clusters.tmp
	mv $CLUSTER/$id.map.clusters.tmp $CLUSTER/$id.map.clusters
	gzip $CLUSTER/$id.map.clusters

	if [ ! -f $MAP/$id.bed.sorted.gz ];
	then 
		gzip $MAP/$id.bed.sorted
	fi

	#### Step-3: flag annotated block groups ####
	if [ -f $CLUSTER/$id.map.clusters.gz ];
	then
		gunzip $CLUSTER/$id.map.clusters.gz
	fi
	if [ ! -f $CLUSTER/$id.map.clusters ];
	then
		echo "Cannot locate $CLUSTER/$id.map.clusters"
		exit		
	fi

	flagCluster.pl -a $ANNOTATION/ncRNAs.hg19.bed.format -c $CLUSTER/$id.map.clusters -f 0 > $CLUSTER/$id.map.clusters.flagged
	flagCluster.pl -a $ANNOTATION/loci.hg19.bed.format -c $CLUSTER/$id.map.clusters.flagged -f 1 -o 50 > $CLUSTER/$id.map.clusters.flagged.tmp
	mv $CLUSTER/$id.map.clusters.flagged.tmp $CLUSTER/$id.map.clusters.flagged

	if [ ! -f $CLUSTER/$id.map.clusters.gz ];
	then
		gzip $CLUSTER/$id.map.clusters
	fi

	#### Step-4: extract significant and annotated block groups ####
	if [ ! -f $CLUSTER/$id.map.clusters.flagged ]
	then
		echo "Cannot locate $CLUSTER/$id.map.clusters.flagged"
		exit
	fi
	echo "Extracting significant block groups for $id"
	if [ -f $CLUSTER/$id.map.clusters.flagged.sig ];
	then
		rm $CLUSTER/$id.map.clusters.flagged.sig
	fi
	perl -an -F'/\s+/' -e 'if($_=~/^>/) { $isSig=0; $size=($F[3]-$F[2])+1; if($size>=50 && $size<=200 && $F[7]>1) { print $_; $isSig=1; } } elsif($isSig) { print $_; }' $CLUSTER/$id.map.clusters.flagged >> $CLUSTER/$id.map.clusters.flagged.sig
	if [ -f $CLUSTER/$id.map.clusters.flagged.sig.anno ];
	then
		rm $CLUSTER/$id.map.clusters.flagged.sig.anno
	fi
	perl -an -F'/\s+/' -e 'if($_=~/^>/) { $isSig=0; $size=($F[3]-$F[2])+1; if($size>=50 && $size<=200 && $F[7]>1 && $F[9]!~/n\/a/) { print $_; $isSig=1; } } elsif($isSig) { print $_; }' $CLUSTER/$id.map.clusters.flagged >> $CLUSTER/$id.map.clusters.flagged.sig.anno

	#### Step-5: compute block group and block statistics
	if [ ! -f $CLUSTER/bgAnnoDist.out ];
	then
		perl -e 'print "#id\tmiRNA\tsnoRNA\ttRNA\tOthers\tUnannotated\tTotal\n";' >> $CLUSTER/bgAnnoDist.out
	fi

	grep "^>" $CLUSTER/$id.map.clusters.flagged | perl -an -F'/\s+/' -e '$total++; if($_=~/miRNA/) { $mirna++ } elsif($_=~/snoRNA/) { $snorna++; } elsif($_=~/tRNA/) { $trna++; } elsif($_=~/n\/a/) { $unannotated++ } else { $others++; } END { print "BG_'$id'_All\t$mirna\t$snorna\t$trna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	grep "^>" $CLUSTER/$id.map.clusters.flagged.sig | perl -an -F'/\s+/' -e '$total++; if($_=~/miRNA/) { $mirna++ } elsif($_=~/snoRNA/) { $snorna++; } elsif($_=~/tRNA/) { $trna++; } elsif($_=~/n\/a/) { $unannotated++ } else { $others++; } END { print "BG_'$id'_Sig\t$mirna\t$snorna\t$trna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	grep "^>" $CLUSTER/$id.map.clusters.flagged.sig.anno | perl -an -F'/\s+/' -e '$total++; if($_=~/miRNA/) { $mirna++ } elsif($_=~/snoRNA/) { $snorna++; } elsif($_=~/tRNA/) { $trna++; } elsif($_=~/n\/a/) { $unannotated++ } else { $others++; } END { print "BG_'$id'_Anno\t$mirna\t$snorna\t$trna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	echo "Done...."

	compBlockStat.pl -c $CLUSTER/$id.map.clusters.flagged -p 3 > $BLOCKSTAT/$id.map.clusters.flagged.bgStat
	compBlockStat.pl -c $CLUSTER/$id.map.clusters.flagged -p 2 > $BLOCKSTAT/$id.map.clusters.flagged.bStat
	R --no-save --vanilla --slave < $SOFTWARE/plotBlockStat.r --args $id $BLOCKSTAT/$id.map.clusters.flagged.bgStat $BLOCKSTAT/$id.map.clusters.flagged.bStat $BLOCKSTAT

	#### Step-6: convert blockbuster to ucsc format ####
	if [ -f $CLUSTER/$id.map.clusters.flagged ];
	then
		blockbuster2ucsc.pl -b $CLUSTER/$id.map.clusters.flagged -l 1 -t $id -o $UCSC/$id.map.clusters.flagged
	fi
done
