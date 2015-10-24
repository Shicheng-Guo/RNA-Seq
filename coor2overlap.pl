#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

##########################################################################################################
# parse input options
use vars qw($coor1 $coor2 $strand1 $strand2 $perOverlap $help);
$perOverlap=1;

GetOptions ("c1=s" => \$coor1,
            "c2=s" => \$coor2,
            "s1=s" => \$strand1,
            "s2=s" => \$strand2,
            "p"    => \$perOverlap,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$coor1 || !$coor2 || !$strand1 || !$strand2);

##########################################################################################################
sub usage {
	print STDERR "\nProgram: coor2overlap.pl (determine overlap between two coordinates)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2overlap.pl -c1 <coordinate> -c2 <coordinate> -s1 <strand> -s2 <strand> [OPTIONS]\n";
	print STDERR " -c1 <coordinate>   [chr:start-end]\n";
	print STDERR " -c2 <coordinate>   [chr:start-end]\n";
	print STDERR " -s1 <strand>       [+/-]\n";
	print STDERR " -s2 <strand>       [+/-]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -p                 [overlap in percentage (default: 1)]\n\n";
	exit(-1);
}
##########################################################################################################

my $overlap=0;

my($chr1, $start1, $end1)=split(/[\:\-]+/, $coor1);
my($chr2, $start2, $end2)=split(/[\:\-]+/, $coor2);

$overlap=checkOverlap($start1, $end1, $start2, $end2, $perOverlap, $chr1, $chr2, $strand1, $strand2);

print "$overlap\n";

exit;
