#!/usr/bin/perl -w
##########################################################################################
## perl program to intersect two bed files
##########################################################################################

use strict;
use warnings;
use Getopt::Long;
use perlModule;

#############################################################################################################
use vars qw($bedFile1 $bedFile2 $help $overlapThreshold);
$overlapThreshold=1;

GetOptions ("f=s"  => \$bedFile1,
            "s=s"  => \$bedFile2,
            "o=i"  => \$overlapThreshold,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$bedFile1 || !$bedFile2);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: bed2anno.pl (retrieve annotations (bed) overlapping to query bed file)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: bed2intersect.pl -q <bed file> -a <bed file>\n";
	print STDERR " -f <query bed file>\n";
	print STDERR " -s <annotation bed file>\n";
	print STDERR " -o <overlap>          [default 1 nt]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}

#############################################################################################################

## open first bed file
open(INFILE, "$bedFile1") || die $!;
my @bedQ=<INFILE>;
close INFILE;

## open second bed file
open(INFILE, "$bedFile2") || die $!;
my @bedA=<INFILE>;
close INFILE;

foreach my $l1(@bedQ) {
	next if($l1=~/^\#/);
	my @q=();
	parseBED($l1, \@q);
	foreach my $l2(@bedA) {
		next if($l2=~/^\#/);
		my @a=();
		parseBED($l2, \@a);
		next if($q[0]!~/^$a[0]$/);
		next if($q[5]!~/^$a[5]$/);
		my $overlap = checkOverlap($q[1], $q[2], $a[1], $a[2], 0);

		if($overlap>=$overlapThreshold) {
			print join "\t",@q; print "\t";
			if(defined($a[7])) { print "$a[7],"; }
			else { print "$a[3],"; }
			#print join "\t",@s; print "\n";
			#last;
		}
	}
	print "\n";
}

exit(0);
