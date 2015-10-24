#!/bin/bash
#PBS -l nodes=1:ppn=4

ORG="hsa"
ORGDIR="humanBodyMap"
#ORGDIR="emtab305"
#ORGDIR="emtab305_trim"
#ORGDIR1="emtab305_norm"
#ORGDIR="egeod28123"

READ=/home/projects/deepalign/data/reads/$ORG/$ORGDIR/raw
TAG=/home/projects/deepalign/data/reads/$ORG/$ORGDIR/tags
MAP=/home/projects/deepalign/data/maps/$ORG/$ORGDIR
CLUSTER=/home/projects/deepalign/data/clusters/$ORG/$ORGDIR
UCSC=/home/projects/deepalign/data/ucsc/$ORG/$ORGDIR
BLOCKSTAT=/home/projects/deepalign/data/compBlockStat/$ORG/$ORGDIR
ANNOTATION=/home/projects/deepalign/data/annotations/$ORG
SOFTWARE=/home/users/sachin/software/myScripts
FASTQ=$1

if [ ! -d "$READ" -o ! -d "$TAG" -o ! -d "$MAP" -o ! -d "$CLUSTER" -o ! -d "$UCSC" -o ! -d "$BLOCKSTAT" -o ! -d "$ANNOTATION" -o ! -d "$SOFTWARE" ];
then
	echo "Cannot locate one or more input directories"
	exit
fi

for fastq in $READ/*.fastq.gz
#for fastq in "${FASTQ[@]}"
do
	id=`echo $fastq | perl -ane '@L=split(/\//,$F[0]);foreach(@L){if(/fastq\.gz/){($id)=split(/\./,$_);print "$id\n";}}'`

comment1()
{
exit
	#### Step-1: filter (Clipped-Trimmed-Quality, -Q 33 specify that read data is in standard SANGER format) ####
	if [ ! -f $TAG/$id.fasta.gz -a "$ORGDIR" = "humanBodyMap" ];
	then
		echo "Performing analysis for $fastq ($id)" > $TAG/$id.log
		zcat $fastq | fastx_clipper -a ATCTCGTATGCCGTCTTCTGCTTG -Q 33 -v | fastx_clipper -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACG -Q 33 -v | fastx_trimmer -Q 33 -f 6 -l 45 | fastx_trimmer -m 18 -Q 33 -v | fastq_quality_trimmer -t 20 -l 18 -Q 33 -v | fastq_quality_filter -q 20 -p 80 -o $TAG/$id.filter -Q 33 -v >> $TAG/$id.log 

		## format (Collapsed-Segemehl compatible Fasta format) ####
		fastq_to_fasta -r -n -i $TAG/$id.filter -Q 33 | grep -v "^>" | sort | uniq -c | perl -an -F'/\s+/' -e '$i++; print ">'$id'_$i|$F[1]\n$F[2]\n";' > $TAG/$id.fasta
		echo "Done, quality filter for $fastq ($id)" >> $TAG/$id.log
		gzip $TAG/$id.fasta
	elif [ ! -f $TAG/$id.fasta.gz -a "$ORGDIR" = "emtab305" ];
	then
		echo "Performing analysis for $fastq ($id)" > $TAG/$id.log
		zcat $fastq | fastx_trimmer -m 18 -v | fastq_quality_trimmer -t 20 -l 18 -v | fastq_quality_filter -q 20 -p 80 -o $TAG/$id.filter -v >> $TAG/$id.log 

		## format (Collapsed-Segemehl compatible Fasta format) ####
		fastq_to_fasta -r -n -i $TAG/$id.filter | grep -v "^>" | sort | uniq -c | perl -an -F'/\s+/' -e '$seq=reverse($F[2]); $seq=~tr/[ATGC]/[TACG]/; $i++; print ">'$id'_$i|$F[1]\n$seq\n";' > $TAG/$id.fasta
		echo "Done, quality filter for $fastq ($id)" >> $TAG/$id.log
		gzip $TAG/$id.fasta
	elif [ ! -f $TAG/$id.fasta.gz -a "$ORGDIR" = "emtab305_trim" ];
	then
		echo "Performing analysis for $fastq ($id)" > $TAG/$id.log
		zcat $fastq | fastx_trimmer -f 8 -l 35 | fastx_trimmer -m 18 -v | fastq_quality_trimmer -t 20 -l 18 -v | fastq_quality_filter -q 20 -p 80 -o $TAG/$id.filter -v >> $TAG/$id.log 

		## format (Collapsed-Segemehl compatible Fasta format) ####
		fastq_to_fasta -r -n -i $TAG/$id.filter | grep -v "^>" | sort | uniq -c | perl -an -F'/\s+/' -e '$seq=reverse($F[2]); $seq=~tr/[ATGC]/[TACG]/; $i++; print ">'$id'_$i|$F[1]\n$seq\n";' > $TAG/$id.fasta
		echo "Done, quality filter for $fastq ($id)" >> $TAG/$id.log
		gzip $TAG/$id.fasta
	elif [ ! -f $TAG/$id.fasta.gz -a "$ORGDIR" = "egeod28123" ];
	then
		echo "Performing analysis for $fastq ($id) for $ORGDIR" > $TAG/$id.log
		zcat $fastq > $TAG/$id.filter
		#zcat $fastq | fastx_trimmer -m 18 -Q 33 -v | fastq_quality_trimmer -t 20 -l 18 -Q 33 -v | fastq_quality_filter -q 20 -p 80 -o $TAG/$id.filter -Q 33 -v >> $TAG/$id.log 

		## format (Collapsed-Segemehl compatible Fasta format) ####
		fastq_to_fasta -r -n -i $TAG/$id.filter -Q 33 | grep -v "^>" | sort | uniq -c | perl -an -F'/\s+/' -e '$i++; print ">'$id'_$i|$F[1]\n$F[2]\n";' > $TAG/$id.fasta
		echo "Done, quality filter for $fastq ($id)" >> $TAG/$id.log
		gzip $TAG/$id.fasta
	fi

	#### Step-2: mapping ####
	if [ ! -f $MAP/$id.map.gz ];
	then
		if [ -f $TAG/$id.fasta.gz ];
		then
			gunzip $TAG/$id.fasta.gz
		fi
		if [ ! -f $TAG/$id.fasta ];
		then
			echo "Cannot locate $TAG/$id.fasta"
			exit		
		fi
		segemehl.x -s --minsize 18 -H 2 -t 8 -K -i /chromo/rna01/nobackup/anthon/hg19/hg19.idx -d /chromo/rna01/nobackup/anthon/hg19/hg19.fa -q $TAG/$id.fasta | gzip > $MAP/$id.map.gz
		if [ ! -f $TAG/$id.fasta.gz ];
		then
			gzip $TAG/$id.fasta
		fi
	fi

	#### Step-3: reformat map file (sort by strand, start, end and chromosome) ####
	if [ ! -f $MAP/$id.map.sorted.gz -a ! -f $MAP/$id.map.sorted ];
	then
		if [ ! -f $MAP/$id.map.gz ];
		then
			echo "Cannot locate $MAP/$id.map.gz"
			exit
		fi
		# uncomment, if using segemehl version prior to 0_1
		#zcat $MAP/$id.map.gz | grep -v "\#" | grep -v -P "^:" | grep -v -P "^\/" | sed 's/>//g' | perl -an -F'/\|/' -e 'chomp($_); if(!defined($seen{$F[0]})) { %seen=(); $seen{$F[0]}=1; foreach $line(@freq) { print "$line\t".scalar(@freq)."\n"; } @freq=(); push(@freq, $_); } else { push(@freq, $_); }  END { foreach $line(@freq) { print "$line\t".scalar(@freq)."\n"; } }'> $MAP/$id.map.sorted
		#cat $MAP/$id.map.sorted | sort -k 11,11 -k 14,14 -k 12n,12 -k 13n,13 -o $MAP/$id.map.tmp

		zcat $MAP/$id.map.gz | grep -v "\#" | grep -v -P "^:" | grep -v -P "^\/" > $MAP/$id.map.sorted
		cat $MAP/$id.map.sorted | sort -k 10,10 -k 13,13 -k 11n,11 -k 12n,12 -o $MAP/$id.map.tmp
		mv $MAP/$id.map.tmp $MAP/$id.map.sorted
		gzip $MAP/$id.map.sorted
		#### (OPTIONAL): convert map to bed format ####
		#zcat $MAP/$id.map.sorted.gz | perl -an -F'/\s+/' -e '($id, $freq)=split(/\|/, $F[0]); $normExp=$freq/$F[15]; print "$F[13]\t$F[11]\t$F[12]\t$id\t$normExp\t$F[10]\n";' > $MAP/$id.map.sorted.bed
	fi
}

	#### Step-4: compute read and tag statistics ####
	perl -e 'print "#id\treads_pre_filter\tread_post_filter\tmapped\tmapped_uniq\n";' >> $TAG/reads.stat
	count=`zcat $fastq | grep "^+" | wc -l`
	perl -e 'print "'$id'\t'$count'\t";' >> $TAG/reads.stat
	count=`zcat $TAG/$id.fasta.gz | grep "^>" | perl -an -F'/\|/' -e '$sum+=$F[1]; END { print "$sum\n"; }'`
	perl -e 'print "'$count'\t";' >> $TAG/reads.stat
	count=`zless $MAP/$id.map.gz | grep -v "\#" | cut -f 2 | sort | uniq | perl -an -F'/\|/' -e '$sum+=$F[1]; END { print "$sum\n"; }'`
	perl -e 'print "'$count'\t";' >> $TAG/reads.stat
	count=`zless $MAP/$id.map.gz | grep -v "\#" | perl -an -F'/\s+/' -e 'if($F[14]==1) { @t=split(/\|/, $F[1]); $sum+=$t[1]; } END { print "$sum\n"; }'`
	perl -e 'print "'$count'\n";' >> $TAG/reads.stat
	cat $TAG/reads.stat | perl -an -F'/\s+/' -e 'chomp($_); if($_=~/^\#/ || $_=~/\(/) { print "$_\n"; } else { $map=sprintf("%0.1f", ($F[3]/$F[2])*100); $uniqMap=sprintf("%0.1f", ($F[4]/$F[3])*100); print "$F[0]\t$F[1]\t$F[2]\t$F[3]($map)\t$F[4]($uniqMap)\n"; }' > $TAG/reads.stat.tmp
	mv $TAG/reads.stat.tmp $TAG/reads.stat

comment2()
{
exit;
	#### Step-5: define block and block groups for each chromosome using blockbuster ####
	#if [ "$id" = "ERR015534Adipose" -o  "$id" = "ERR015537Hypothalamus" -o "$id" = "ERR015539Liver" ];
	#then
	#	blockDist=50
	#else
	#	blockDist=70
	#fi
	#blockDist=40
	blockDist=50

	if [ -f $MAP/$id.map.sorted.gz ];
	then
		gunzip $MAP/$id.map.sorted.gz
	fi
	if [ ! -f $MAP/$id.map.sorted ];
	then
		echo "Cannot locate $MAP/$id.map.sorted"
		exit		
	fi
	CHR=(`cut -f 13 $MAP/$id.map.sorted | grep -v "_" | uniq | sort | uniq`)

	## remove output file, if exists
	if [ -f $CLUSTER/$id.map.clusters.gz ];
	then
		rm $CLUSTER/$id.map.clusters.gz
	fi
	if [ -f $CLUSTER/$id.map.clusters ];
	then
		rm $CLUSTER/$id.map.clusters
	fi

	for chr in "${CHR[@]}"
	do
		#blockbuster.x -format 2 -minBlockHeight 10 -minClusterHeight 50 -scale 0.5 -chr $chr -print 2 $MAP/$id.map.sorted >> $CLUSTER/$id.map.clusters
		#blockbuster.x -format 2 -distance $blockDist -minBlockHeight 2 -minClusterHeight 10 -scale 0.5 -chr $chr -print 2 $MAP/$id.map.sorted >> $CLUSTER/$id.map.clusters
		blockbuster.x -format 2 -distance $blockDist -minBlockHeight 10 -minClusterHeight 10 -scale 0.5 -chr $chr -print 2 -blockHeight rel $MAP/$id.map.sorted >> $CLUSTER/$id.map.clusters
	done

	## validate block groups defined by blockbuster
	## IMPORTANT NOTE: useful for correct ordering of block groups
	validateBlockbuster.pl -c $CLUSTER/$id.map.clusters -d 1000000 -o 100 > $CLUSTER/$id.map.clusters.tmp
	mv $CLUSTER/$id.map.clusters.tmp $CLUSTER/$id.map.clusters
	gzip $CLUSTER/$id.map.clusters

	if [ ! -f $MAP/$id.map.sorted.gz ];
	then 
		gzip $MAP/$id.map.sorted
	fi

	#### Step-6: flag annotated block groups ####
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
	flagCluster.pl -a $ANNOTATION/rfam10.1.hg19.bed.format -c $CLUSTER/$id.map.clusters.flagged -f 0 > $CLUSTER/$id.map.clusters.flagged.tmp
	mv $CLUSTER/$id.map.clusters.flagged.tmp $CLUSTER/$id.map.clusters.flagged
	flagCluster.pl -a $ANNOTATION/loci.hg19.bed.format -c $CLUSTER/$id.map.clusters.flagged -f 1 -o 50 > $CLUSTER/$id.map.clusters.flagged.tmp
	mv $CLUSTER/$id.map.clusters.flagged.tmp $CLUSTER/$id.map.clusters.flagged

	if [ ! -f $CLUSTER/$id.map.clusters.gz ];
	then
		gzip $CLUSTER/$id.map.clusters
	fi

	#### Step-7: extract significant and annotated block groups ####
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
	perl -an -F'/\s+/' -e 'if($_=~/^>/) { $isSig=0; $size=($F[3]-$F[2])+1; if($size<=500) { print $_; $isSig=1; } } elsif($isSig) { print $_; }' $CLUSTER/$id.map.clusters.flagged >> $CLUSTER/$id.map.clusters.flagged.sig

	if [ -f $CLUSTER/$id.map.clusters.flagged.sig.anno ];
	then
		rm $CLUSTER/$id.map.clusters.flagged.sig.anno
	fi
	perl -an -F'/\s+/' -e 'if($_=~/^>/) { $isSig=0; $size=($F[3]-$F[2])+1; if($size<=500 && $F[9]!~/n\/a/) { print $_; $isSig=1; } } elsif($isSig) { print $_; }' $CLUSTER/$id.map.clusters.flagged >> $CLUSTER/$id.map.clusters.flagged.sig.anno

	#### Step-8: compute block group and block statistics
	perl -e 'print "#id\tmiRNA\tsnoRNA\ttRNA\tsnRNA\tscRNA\tyRNA\tOthers\tUnannotated\tTotal\n";' >> $CLUSTER/bgAnnoDist.out
	grep "^>" $CLUSTER/$id.map.clusters.flagged | perl -an -F'/\s+/' -e '$total++; if($F[9]=~/miRNA/) { $mirna++ } elsif($F[9]=~/snoRNA/) { $snorna++; } elsif($F[9]=~/tRNA/) { $trna++; } elsif($F[9]=~/snRNA/) { $snrna++; } elsif($F[9]=~/scRNA/) { $scrna++; } elsif($F[9]=~/Y_RNA/) { $yrna++; } elsif($F[9]=~/n\/a/) { $unannotated++; } else { $others++; } END { print "BG_'$id'_all\t$mirna\t$snorna\t$trna\t$snrna\t$scrna\t$yrna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	grep "^>" $CLUSTER/$id.map.clusters.flagged.sig | perl -an -F'/\s+/' -e '$total++; if($F[9]=~/miRNA/) { $mirna++ } elsif($F[9]=~/snoRNA/) { $snorna++; } elsif($F[9]=~/tRNA/) { $trna++; } elsif($F[9]=~/snRNA/) { $snrna++; } elsif($F[9]=~/scRNA/) { $scrna++; } elsif($F[9]=~/Y_RNA/) { $yrna++; } elsif($F[9]=~/n\/a/) { $unannotated++; } else { $others++; } END { print "BG_'$id'_sig\t$mirna\t$snorna\t$trna\t$snrna\t$scrna\t$yrna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	grep "^>" $CLUSTER/$id.map.clusters.flagged.sig.anno | perl -an -F'/\s+/' -e '$total++; if($F[9]=~/miRNA/) { $mirna++ } elsif($F[9]=~/snoRNA/) { $snorna++; } elsif($F[9]=~/tRNA/) { $trna++; } elsif($F[9]=~/snRNA/) { $snrna++; } elsif($F[9]=~/scRNA/) { $scrna++; } elsif($F[9]=~/Y_RNA/) { $yrna++; } elsif($F[9]=~/n\/a/) { $unannotated++; } else { $others++; } END { print "BG_'$id'_anno\t$mirna\t$snorna\t$trna\t$snrna\t$scrna\t$yrna\t$others\t$unannotated\t$total\n"; }' >> $CLUSTER/bgAnnoDist.out
	echo "Done...."

	compBlockStat.pl -c $CLUSTER/$id.map.clusters.flagged -p 3 > $BLOCKSTAT/$id.map.clusters.flagged.bgStat
	compBlockStat.pl -c $CLUSTER/$id.map.clusters.flagged -p 2 > $BLOCKSTAT/$id.map.clusters.flagged.bStat

	R --no-save --vanilla --slave < $SOFTWARE/plotBlockStat.r --args $id $BLOCKSTAT/$id.map.clusters.flagged.bgStat $BLOCKSTAT/$id.map.clusters.flagged.bStat $BLOCKSTAT
}

comment3()
{
exit
	#### Step-9: convert blockbuster to ucsc format ####
	if [ -f $CLUSTER/$id.map.clusters.flagged ];
	then
		blockbuster2ucsc.pl -b $CLUSTER/$id.map.clusters.flagged -l 1 -t $id -o $UCSC/$id.map.clusters.flagged
	fi
}
done
