#!/bin/sh

usage()
{
  echo >&2 "Usage: $0 -i idxdir -q query [-w workdir] dbfiles"
  echo >&2 "Match a query file with a set of segemehl indices"
  echo >&2 "[OPTIONS]"
  echo >&2 "-D  number of differences"
  echo >&2 "-T  number of threads"
  echo >&2 "-A  accuracy"
  echo >&2 "-B  best only"
  exit 1
}

segemehldir=
segemehltools=
reporttype=
workdir= 
wflag=
idxdflag=
idxdir=
gflag=
mflag=
oflag=
qflag=
xflag=
queryfile=
optstr=
aflag=
diff=1
annotationfilen=
outfiles=
gluefile=
glueopt=

qbname=
acc=""
bestonly=""
thd=""

if [ -z $segemehltools ]; then
  segemehltools="./tools/glue"
fi

while getopts si:a:gR:o:w:D:q:T:A:E: o
do  case "$o" in
  i)    idxdflag=1
  idxdir="$OPTARG";;
  w)    wflag=1   
  workdir="$OPTARG";;
  q)    qflag=1    
  queryfile="$OPTARG";;
  D)    diff="$OPTARG";;
  A)    acc="-A $OPTARG";;
  T)    thd="--threads $OPTARG";;
  B)    bestonly="--bestonly";;

  [?])  usage;;
esac
done
shift $(($OPTIND-1))

if [ ! -z "$idxdflag" ]; then
  printf 'idx directory %s specified.\n' "$idxdir"
fi

if [ $# -eq 0 ]; then
  usage
fi

if [ -z "$qflag" ]; then
  echo >&2 "no query specified."i
  usage
fi



suffix=`echo $1 | grep -o \.fa[sta]*$` 
qbname=`basename $queryfile $suffix`

while [ $# -gt 0 ] 
do
  bname=
  out=
  suffix=`echo $1 | grep -o \.fa[sta]*$` 
  bname=`basename $1 $suffix`


  if [ -z "$bname" -o -z "$suffix" ]; then
    echo >&2 "$1 doesn't seem to have a proper extension (.fa[sta]*)"
    echo >&2 "Exit forced."
    exit 1
  fi

  if [ ! -z "$idxdflag" ]; then
    idxname="$idxdir/$bname.idx" 
  else
    idxname="$bname.idx"
  fi

  time ./segemehl.x  -d $1 -i $idxname -q $queryfile -D $diff $acc $thd > $workdir/$bname.$qbname.mat  

  shift
done

#end of the show 
