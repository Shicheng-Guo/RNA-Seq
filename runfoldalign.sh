#!/bin/bash
#PBS -l nodes=1:ppn=4

SHUFFLE=10
COMPUTE=0
THRESHOLD=300
HELP=0

#### usage ####
usage() {
	echo
	echo Program: "runfoldalign.sh (run foldalign and evaluate the alignments)"
	echo Author: RTH, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: sachin@rth.dk
	echo "Usage: runfoldalign.sh -f <coor1> -s <coor2> -l <shufflings> -c <compute> -t <threshold>"
	echo "Options:"
	echo " -f <coor1>        coordinate of first block group (chr:start-end#strand)"
	echo " -s <coor1>        coordinate of second block group (chr:start-end#strand)"
	echo "Optional:"
	echo " -l <shufflings>   default (10)"
	echo " -c <compute>      [0: alignment, 1: p-value] default(0)"
	echo " -t <threshold>    default(300)"
	echo " -h                help"
	echo
	exit 0
}

#### parse options
while getopts f:s:l:c:t:h: ARG; do
	case "$ARG" in
		f) COOR1=$OPTARG;;
		s) COOR2=$OPTARG;;
		l) SHUFFLE=$OPTARG;;
		c) COMPUTE=$OPTARG;;
		t) THRESHOLD=$OPTARG;;
		h) HELP=1;;
	esac
done

if [ ! $COOR1 -o ! $COOR2 -o $HELP -eq 1 ]; then
	usage
fi

echo "# running foldalign for COOR1: $COOR1, COOR2: $COOR2"
echo "# arguments SHUFFLE: $SHUFFLE, COMPUTE: $COMPUTE, THRESHOLD: $THRESHOLD"

STATUS=`echo "$COOR1 $COOR2" | perl -an -F'/\s+/' -e '($coor1, $strand1)=split(/\#/,$F[0]); ($coor2, $strand2)=split(/\#/,$F[1]); if(!defined($coor1) || !defined($coor2) || !defined($strand1) || !defined($strand2)) { print -1; } else { print 0; }'`

if [ $STATUS -eq -1 ]; then
	usage
fi

t1=`echo $COOR1 | perl -ane '$_=~s/\#.+//g; $_=~s/[\:\-]+//g; print $_;'`
t2=`echo $COOR2 | perl -ane '$_=~s/\#.+//g; $_=~s/[\:\-]+//g; print $_;'`

FILE1=`echo ${t1}_${t2}_1`
FILE2=`echo ${t1}_${t2}_2`

if [ -e $FILE1 ]; then
	rm $FILE1;
fi
if [ -e $FILE2 ]; then
	rm $FILE2;
fi

if [ $COMPUTE -eq 0 ]; then
	echo "$COOR1 $COOR2 $FILE1 $FILE2" | perl -an -F'/\s+/' -e '($coor1, $strand1)=split(/\#/,$F[0]); ($coor2, $strand2)=split(/\#/,$F[1]); system("coor2seq.pl -c $coor1 -s $strand1 -t unannotated > $F[2]"); system("coor2seq.pl -c $coor2 -s $strand2 -t annotated > $F[3]"); system("foldalign $F[2] $F[3] | grep -w ALIGN");';
	rm $FILE1
	rm $FILE2

elif [ $COMPUTE -eq 1 ]; then
	ACTSCORE=`echo "$COOR1 $COOR2 $FILE1 $FILE2" | perl -an -F'/\s+/' -e '($coor1, $strand1)=split(/\#/,$F[0]); ($coor2, $strand2)=split(/\#/,$F[1]); system("coor2seq.pl -c $coor1 -s $strand1 > $F[2]"); system("coor2seq.pl -c $coor2 -s $strand2 > $F[3]"); @foldalign=\`foldalign $F[2] $F[3] | grep -w "ALIGN"\`; foreach(@foldalign) { if($_=~/Score\:/) { $_=~s/^.+\:\s+//g; print $_; } }'`
	if [ $ACTSCORE -gt $THRESHOLD ]; then
		shuffle -n $SHUFFLE -d $FILE1 | perl -ane 'BEGIN { $i=1; } if($_=~/^>/) { $_=~s/^\>[^\s]+/\>shuffle_seq$i/g; print $_; $i++; } else { print $_; }' > $FILE1.shuffle
		shuffle -n $SHUFFLE -d $FILE2 | perl -ane 'BEGIN { $i=1; } if($_=~/^>/) { $_=~s/^\>[^\s]+/\>shuffle_seq$i/g; print $_; $i++; } else { print $_; }' > $FILE2.shuffle
		FP=0;
		for ((i=1; i<=SHUFFLE; i++)); do
			for ((j=1; j<=SHUFFLE; j++)); do
				faOneRecord $FILE1.shuffle shuffle_seq$i > $FILE1.shuffle.tmp;
				faOneRecord $FILE2.shuffle shuffle_seq$j > $FILE2.shuffle.tmp;
				SHUFFLESCORE=`foldalign $FILE1.shuffle.tmp $FILE2.shuffle.tmp | grep -w "Score:" | sed 's/ //g' | sed 's/;ALIGNScore://'`
				#echo "Shuff. score: $SHUFFLESCORE"
				if [ $SHUFFLESCORE -gt $ACTSCORE ]; then
					FP=$((FP+1))
				fi
			done
		done
		FP=$((FP+1))

		# empirical p-value = (r+1)/(n+1) [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/]
		PVALUE=`echo "$FP $SHUFFLE" | perl -an -F'/\s+/' -e 'printf("%0.2f", $F[0]/(($F[1]*$F[1])+1));'`
		echo "$COOR1 $COOR2 $ACTSCORE $PVALUE"
		rm $FILE1.shuffle
		rm $FILE2.shuffle
		rm $FILE1.shuffle.tmp
		rm $FILE2.shuffle.tmp
	else
		echo "$COOR1 $COOR2 $ACTSCORE NS"
	fi
	rm $FILE1
	rm $FILE2
else
	usage
fi
