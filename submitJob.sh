#!/bin/bash
#PBS -l nodes=1:ppn=4

## EDIT
OPTION=1

## EDIT
#ORGDIR="sperm"
ORGDIR="shortTotalRNA"
#ORGDIR="emtab305_norm"
#ORGDIR="egeod28123"

## EDIT
OUTDIR="/home/projects/deepalign/data/tree/hsa/$ORGDIR/UNIQ/all/"
PROGDIR="~/software/myScripts"

if [ $OPTION -eq 0 ]
then
	## EDIT
	INDIR="/home/projects/deepalign/data/expAna/hsa/$ORGDIR"
	echo "$PROGDIR/rnaExpAna.pl -c a -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR -f $INDIR/bgExpCons" | qsub -q default
elif [ $OPTION -eq 1 ]
then
	## EDIT
	INDIR="/home/projects/deepalign/data/expAna/hsa/$ORGDIR/UNIQ/statTmp"
	for infile in $INDIR/*Tmp*
	do
		## shortTotalRNA
		#echo "$PROGDIR/rnaExpAnaAll.pl -c b -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR/ -f $infile -o $OUTDIR -r /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/reads.stat -m /home/projects/deepalign/data/maps/hsa/$ORGDIR/index/ -z /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/sizeFactorInd -q > $infile.stat" | qsub -q fastlane
		echo "$PROGDIR/rnaExpAnaAll.pl -c b -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR/ -f $infile -o $OUTDIR -r /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/reads.stat -m /home/projects/deepalign/data/maps/hsa/$ORGDIR/index/ -z /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/sizeFactorInd -e -q > $infile.stat" | qsub -q fastlane
		## egeod28123
		#echo "$PROGDIR/rnaExpAnaAll.pl -c b -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR/ -f $infile -o $OUTDIR -r /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/reads.stat -m /home/projects/deepalign/data/maps/hsa/$ORGDIR/index/ -z /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/sizeFactor1 -q > $infile.stat" | qsub -q fastlane
		## emtab305_norm
		#echo "$PROGDIR/rnaExpAnaAll.pl -c b -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR/ -f $infile -o $OUTDIR -r /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/reads.stat -m /home/projects/deepalign/data/maps/hsa/emtab305_trim/index/ -z /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/sizeFactor -q > $infile.stat" | qsub -q fastlane
		#echo "$PROGDIR/rnaExpAna.pl -c b -i /home/projects/deepalign/data/clusters/hsa/$ORGDIR/ -f $infile -o $OUTDIR -r /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/reads.stat -m /home/projects/deepalign/data/maps/hsa/emtab305_trim/index/ -z /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/sizeFactor > $infile.stat" | qsub -q fastlane
	done
elif [ $OPTION -eq 2 ]
then
	INDIR="/home/projects/deepalign/data/expAna/hsa/$ORGDIR/sigTmpAll"
	for infile in $INDIR/*Tmp*
	do
		## shortTotalRNA
		echo "$PROGDIR/rnaExpAnaAll.pl -c c -s $infile -o $OUTDIR -m /home/projects/deepalign/data/maps/hsa/$ORGDIR/index/ -r /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/reads.stat -z /home/projects/deepalign/data/reads/hsa/$ORGDIR/tags/sizeFactorInd > $infile.sig" | qsub -q fastlane
		## emtab305_norm
		#echo "$PROGDIR/rnaExpAna.pl -c c -s $infile -o $OUTDIR -m /home/projects/deepalign/data/maps/hsa/emtab305_trim/index/ -r /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/reads.stat -z /home/projects/deepalign/data/reads/hsa/emtab305_trim/tags/sizeFactor > $infile.sig" | qsub -q fastlane
	done
fi
