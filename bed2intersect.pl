#!/usr/bin/perl -w
##########################################################################################
## perl program to intersect two bed files
##########################################################################################

use strict;
use warnings;
use Getopt::Long;
use perlModule;

#############################################################################################################
use vars qw($bedFile1 $bedFile2 $coor $strand $overlapThreshold $printNonOverlap $help);
$overlapThreshold=1;

GetOptions ("f=s"  => \$bedFile1,
            "c=s"  => \$coor,
            "t=s"  => \$strand,
            "s=s"  => \$bedFile2,
            "o=i"  => \$overlapThreshold,
            "v"    => \$printNonOverlap,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || ((!$coor && !$strand) && !$bedFile1) || !$bedFile2);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: bed2intersect.pl (intersect two bed files)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: bed2intersect.pl [-f <file> | -c <coor>] -s <file> [OPTIONS]\n";
	print STDERR " -f <file>        [first BED file (optional)]\n";
	print STDERR " -c <coor>        [coordinate (chr:start-end) (optional)]\n";
	print STDERR " -s <file>        [second BED file (required)]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <overlap>     [default 1 nt (optional)]\n";
	print STDERR " -t <strand>      [strand (+/-) (optional)]\n";
	print STDERR " -v               [print only non-overlapping coordinates]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}

#############################################################################################################

## open first bed file or reformat input coordinate
my @bedF=();
if(defined($bedFile1)) {
	open(INFILE, "$bedFile1") || die $!;
	@bedF=<INFILE>;
	close INFILE;
}
elsif(defined($coor)) {
	my @t=split(/[\:\-]/,$coor);
	if(defined($strand) && ($strand=~/\+/ || $strand=~/\-/)) {
		push(@bedF, "$t[0]\t$t[1]\t$t[2]\tNA\t1\t$strand");
	}
	else {
		push(@bedF, "$t[0]\t$t[1]\t$t[2]");
	}
}
else { usage(); }

## open second bed file
open(INFILE, "$bedFile2") || die $!;
my @bedS=<INFILE>;
close INFILE;

my $found=();
foreach my $l1(@bedF) {
	$found=0;
	next if($l1=~/^\#/);
	my @q=();	parseBED($l1, \@q);
	foreach my $l2(@bedS) {
		next if($l2=~/^\#/);
		my @s=();	parseBED($l2, \@s);
		next if($q[0]!~/^$s[0]$/);

		## case 1: BED format is "chr start end"
		if(scalar(@q)==3 || scalar(@s)==3) {
			my $overlap = checkOverlap($q[1], $q[2], $s[1], $s[2], 0, $q[0], $s[0]);
			if($overlap>=$overlapThreshold) {
				if(!defined($printNonOverlap)) {
					print join "\t",@q; print "\t";
					print join "\t",@s; print "\n";
				}
				$found=1;
				last;
			}
		}
		## case 2: BED format is "chr start end name score strand"
		elsif(scalar(@q)>=5 && scalar(@s)>=5) {
			my $overlap = checkOverlap($q[1], $q[2], $s[1], $s[2], 0, $q[0], $s[0], $q[5], $s[5]);
			if($overlap>=$overlapThreshold) {
				if(!defined($printNonOverlap)) {
					print join "\t",@q; print "\t";
					print join "\t",@s; print "\n";
				}
				$found=1;
				last;
			}
		}
		else {
			print STDERR "Incorrect bed format\n\t$l1\t$l2\n"; exit;
		}
	}
	if($found==0) {
		if(defined($printNonOverlap)) {
			print join "\t",@q; print "\n";
		}
		else {
			print "NA\n";
		}
	}
}

exit(0);

