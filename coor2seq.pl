#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use LWP::Simple;
use perlModule;

##########################################################################################################
# parse input options
my %coor=();
use vars qw($help $title) ;
$coor{'strand'}="+";
$coor{'db'}="hg19";
$coor{'source'}="fasta";

GetOptions ("c=s"  => \$coor{'coor'},
            "s=s"  => \$coor{'strand'},
            "d=s"  => \$coor{'db'},
            "r=s"  => \$coor{'source'},
            "t=s"  => \$title,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$coor{'coor'});

##########################################################################################################
sub usage {
	print STDERR "\nProgram: coor2seq.pl (given the coordinate, retrieve nucleotide sequence from UCSC)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2seq.pl -c <coordinate> [OPTIONS]\n";
	print STDERR " -c <coordinate>  [chr:start-end]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -s <strand>      [+/- (default: +)]\n";
	print STDERR " -d <string>      [genome (default: hg19)]\n";
	print STDERR " -r <string>      [fasta or ucsc (default: local)]\n"; 
	print STDERR " -t <string>      [title]\n\n";
	exit(-1);
}
##########################################################################################################

($coor{'chr'}, $coor{'start'}, $coor{'end'}) = split(/[\:\-]+/, $coor{'coor'});

if(!$coor{'chr'} || !$coor{'start'} || !$coor{'end'}) { print STDERR "Incorrect coordinate\n"; usage(); }

if($coor{'strand'}=~/\+/) {
	if($coor{'source'}=~/fasta/) {
		$coor{'url'}="twoBitToFa /chromo/ucsc/gbdb/$coor{'db'}/$coor{'db'}.2bit:$coor{'chr'}:".($coor{'start'}-1)."-$coor{'end'} stdout";
		$coor{'seq'}=`$coor{'url'} | grep -v "^\>"`;
		$coor{'seq'}=">$coor{'chr'}:$coor{'start'}-$coor{'end'}\n".$coor{'seq'};
	}
	else {
		$coor{'url'} = "http://ucsc.genome.ku.dk/cgi-bin/hgc?hgsid=1315&g=htcGetDna2&table=&i=mixed&o=58514087&l=58514087&r=58514192&getDnaPos=$coor{'chr'}%3A$coor{'start'}-$coor{'end'}&db=$coor{'db'}&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA";
		$coor{'seq'} = get($coor{'url'});

		$coor{'seq'}=~s/\<[\/]{0,1}PRE\>//g;
		$coor{'seq'}=~s/^\s*\n+//mg;
	}
}
elsif($coor{'strand'}=~/\-/) {
	if($coor{'source'}=~/fasta/) {
		$coor{'url'}="twoBitToFa /chromo/ucsc/gbdb/$coor{'db'}/$coor{'db'}.2bit:$coor{'chr'}:".($coor{'start'}-1)."-$coor{'end'} stdout";
		$coor{'seq'}=`$coor{'url'} | grep -v "^>"`;
		$coor{'seq'}=reverse_complement($coor{'seq'});
		$coor{'seq'}=">$coor{'chr'}:$coor{'start'}-$coor{'end'}\n".$coor{'seq'};
	}
	else {
		$coor{'url'} = "http://ucsc.genome.ku.dk/cgi-bin/hgc?hgsid=1315&g=htcGetDna2&table=&i=mixed&o=58514087&l=58514087&r=58514192&getDnaPos=$coor{'chr'}%3A$coor{'start'}-$coor{'end'}&db=$coor{'db'}&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&hgSeq.revComp=on&boolshad.hgSeq.revComp=0&submit=get+DNA";
		$coor{'seq'} = get($coor{'url'});

		$coor{'seq'}=~s/\<[\/]{0,1}PRE\>//g;
		$coor{'seq'}=~s/^\s*\n+//mg;
	}
}
else { print STDERR "Incorrect strand format $coor{'strand'}\n"; usage(); }

if(defined($title)) { $coor{'seq'}=~s/\>[^\s]+/\>$title/g; }
print "$coor{'seq'}";

exit;
