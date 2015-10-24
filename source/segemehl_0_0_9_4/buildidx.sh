#!/bin/sh

usage()
{
  echo >&2 "Usage: $0 [-i idxdir] -dbfiles"
  echo >&2 "build segemehl indices for a set of (multiple)fasta files"
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
annotationfilen=
outfiles=
gluefile=
glueopt=


if [ -z $segemehltools ]; then
  segemehltools="./tools/glue"
fi

while getopts si:a:gR:o:w:D:q:t:A:E: o
do  case "$o" in
  i)    idxdflag=1
  idxdir="$OPTARG";;
  [?])  usage;;
esac
done
shift $(($OPTIND-1))

if [ ! -z "$idxdflag" ]; then
  printf 'Working directory %s specified.\n' "$idxdir"
fi

if [ $# -eq 0 ]; then
  usage
fi

while [ $# -gt 0 ] 
do
  bname=
  out=

  suffix=`echo $1 | grep -o \.fa[sta]*$` 
  bname=`basename $1 $suffix`

  if [ -z "$bname" -o -z "$suffix" ]; then
    echo >&2 "$1 doesn't seem to have a proper extension (.fa[sta]*)"
    echo >&2 "Exit forced."
  fi

  if [ ! -z "$idxdflag" ]; then
    idxname="$idxdir/$bname.idx" 
  else
    idxname="$bname.idx"
  fi

  time ./segemehl.x -d $1 -x $idxname   

  shift
done

#end of the show 
