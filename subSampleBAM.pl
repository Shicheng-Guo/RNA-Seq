#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Tie::IxHash;

############################################################################################################
# parse input options
use vars qw($bamFile $coorFile $outDir $help);
$outDir=".";

GetOptions ("i=s"  => \$bamFile,
            "j=s"  => \$coorFile,
						"o=s"  => \$outDir,
            "help" => \$help,
            "h"    => \$help);

usage() if(!$bamFile || !$coorFile || $help);

############################################################################################################
sub usage {
	print STDERR "\nProgram: subSampleBAM.pl (extract reads corresponding to input coordinates from a BAM file)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: subSampleBAM.pl -i <file> -j <file> [OPTIONS]\n";
	print STDERR " -i <file>        [input BAM file]\n";
	print STDERR " -j <file>        [coordinate file in BED format]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <dir>         [output directory to store the sub-sampled BAM file (default: current)]\n";
	print STDERR " -h               [help]\n";
	print STDERR "[REQUIREMENTS]\n";
	print STDERR " a) input BAM file is sorted\n";
	print STDERR "    If not, BAM file can be sorted using samtools (http://samtools.sourceforge.net/samtools.shtml)\n";
	print STDERR "    command to sort BAM file is 'samtools sort aln.bam aln.sorted'\n";
	print STDERR " b) BEDTools is installed on the computer\n";
	print STDERR "    If not, it can be downloaded from http://code.google.com/p/bedtools/\n\n";
	exit(-1);
}
############################################################################################################

if(! -e $bamFile) {
	print "\nError: cannot find $bamFile !!!\n";
	usage();
}

if(! -e $coorFile) {
	print "\nError: cannot find $coorFile !!!\n";
	usage();
}

my @reads=`zless $bamFile | bedtools bamtobed -i stdin | intersectBed -wb -sorted -a $coorFile -b stdin`;
my $inputFields=`zless $coorFile | perl -ane 'print scalar(\@F)."\n"; exit;'`;

foreach(@reads) {
	my @t=split(/\s+/,$_);
	for(my $i=$inputFields; $i<scalar(@t); $i++) {
		print "$t[$i]\t";
	} print "\n";
}

exit;

