#!/bin/bash
#PBS -l nodes=1:ppn=4

#### usage ####
usage() {
	echo Program: "indexMap.sh (create indexing of segemehl or bed map file)"
	echo Author: RTH, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: sachin@rth.dk
	echo "Usage: indexMap.sh -i <file> -f <string> [OPTIONS]"
	echo " -i <file>     [input map or anno file]"
	echo " -f <string>   [format of input file (segemehl, bed and anno)]"
	echo "[OPTIONS]:"
	echo " -o <dir>      [output directory with index files]"
	echo " -b            [use UCSC utility \"bedSplitOnChrom\"]"
	echo
	exit 0
}

#### parse options ####
while getopts i:o:f:b ARG; do
	case "$ARG" in
		i) INFILE=$OPTARG;;
		o) OUTDIR=$OPTARG;;
		f) FORMAT=$OPTARG;;
		b) UTILITY=1;;
	esac
done

if [ ! "$INFILE" -o ! "$FORMAT" ]; then
	usage
fi

## create output directory, if do not exist
if [ -z "$OUTDIR" ]; then
	OUTDIR=`echo $INFILE | perl -ane '$_=~s/^.*\///g; $_=~s/\..+$//g; print $_;'`;
	if [ ! -d "$OUTDIR" ]; then
		mkdir $OUTDIR;
	fi
fi

if [ "$FORMAT" == "bed" -a "$UTILITY" == 1 ]; then
	`bedSplitOnChrom $INFILE $OUTDIR -strand`
else
	if [ "$FORMAT" == "segemehl" ];
	then
		CHR=(`zless $INFILE | cut -f 13 | grep -v "_" | sort | uniq`)
	elif [ "$FORMAT" == "bed" -o "$FORMAT" == "anno" ];
	then
		CHR=(`zless $INFILE | cut -f 1 | grep -v "_" | sort | uniq`)
	else
		echo "unknown format of map file"
		exit 0
	fi

	for chr in "${CHR[@]}";
	do
		fileFormat=`echo $INFILE | perl -ane '\$_=~s/\\.gz//g; \$_=~s/.+\///g; print \$_;'`
		echo "$INFILE $chr $fileFormat";
		if [ $FORMAT == "segemehl" ]; then
			`zgrep -E "$chr\\s+.*\\s+\\+\\s*" $INFILE | gzip > $OUTDIR/$fileFormat.$chr.P.gz`;
			`zgrep -E "$chr\\s+.*\\s+\\-\\s*" $INFILE | gzip > $OUTDIR/$fileFormat.$chr.N.gz`;
		elif [ $FORMAT == "bed" ]; then
			`zgrep -E "$chr\\s+.*\\s+\\+\\s*" $INFILE | perl -an -F'/\\s+/' -e '@t=split(/\_/, $F[3]); $expr=sprintf("%0.0f", $t[scalar(@t)-1]*$F[4]); print "Q\t$t[0]_$t[1]|$expr\t1\t1\t1\t1\t1\t1\t1\t$F[5]\t$F[1]\t$F[2]\t$F[0]\tM1;\t$t[2]\n";' | gzip > $OUTDIR/$fileFormat.$chr.P.gz`;
			`zgrep -E "$chr\\s+.*\\s+\\-\\s*" $INFILE | perl -an -F'/\\s+/' -e '@t=split(/\_/, $F[3]); $expr=sprintf("%0.0f", $t[scalar(@t)-1]*$F[4]); print "Q\t$t[0]_$t[1]|$expr\t1\t1\t1\t1\t1\t1\t1\t$F[5]\t$F[1]\t$F[2]\t$F[0]\tM1;\t$t[2]\n";' | gzip > $OUTDIR/$fileFormat.$chr.N.gz`;
		elif [ $FORMAT == "anno" ]; then
			`zgrep -E "$chr\\s+.*\\s+\\+\\s*" $INFILE | gzip > $OUTDIR/$fileFormat.$chr.P.gz`;
			`zgrep -E "$chr\\s+.*\\s+\\-\\s*" $INFILE | gzip > $OUTDIR/$fileFormat.$chr.N.gz`;
		fi
	done
fi
