#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use perlModule;

##########################################################################################################
## parse input options
use vars qw($bgDir $countFile $mapDir $uniqReadsFilter $chrFilter $fileFilter $compSizeFactor $help);
$uniqReadsFilter=0;
$chrFilter=".+";
$fileFilter=".+";

GetOptions ("d=s"  => \$bgDir,
            "r=s"  => \$countFile,
            "m=s"  => \$mapDir,
            "u=i"  => \$uniqReadsFilter,
            "c=s"  => \$chrFilter,
            "f=s"  => \$fileFilter,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$bgDir || !$countFile);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: computeDESeqMatrix.pl (create DESeq matrix with read counts for each block group)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: computeDESeq.pl -b <dir> -r <file> [OPTIONS]\n";
	print STDERR " -d <dir>          [directory with block groups in blockbuster format]\n";
	print STDERR " -r <file>         [output file with read count at each coor (block group) across all samples]\n";
	print STDERR " -m <dir>          [directory with indexed map files]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -u <int>          [compute read count for block groups with at least given number of unique reads]\n";
	print STDERR " -c <string>       [compute read count for block groups of given chromosome only]\n";
	print STDERR " -f <string>       [consider only those block group files that match the given name]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

##########################################################################################################

opendir(INDIR, $bgDir) || die $!;
my @files=();

## read block group information (coordinate, identifier, chr, strand) from each block group file
foreach my $bgFile(sort { $a cmp $b } readdir(INDIR)) {
	if($bgFile=~/.*$fileFilter.*flagged\.sig.filter.uniq.gz$/) {
	#if(($bgFile=~/.*Blood.*flagged\.sig$/ || $bgFile=~/.*Liver.*flagged\.sig$/) && $bgFile=~/Rep1/) {
		push(@files, $bgFile);
	}
}
close INDIR;

## compute read count for each coordinate (block group) across all samples
compReadCount($bgDir, \@files, $countFile, 0, $uniqReadsFilter, $chrFilter, $mapDir);

exit(0);
