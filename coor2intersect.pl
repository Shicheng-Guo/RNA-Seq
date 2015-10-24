#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

############################################################################################################
# parse input options
my %coor=();
use vars qw($help) ;

GetOptions ("f=s"  => \$coor{'coorfile1'}{'name'},
            "s=s"  => \$coor{'coorfile2'}{'name'},
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$coor{'coorfile1'}{'name'} || !$coor{'coorfile2'}{'name'});

############################################################################################################
sub usage {
	print STDERR "\nProgram: coor2intersect.pl (given two sets of coordinates, retrieve common coordinates)\n";
	print STDERR "Note: output will be from coorfile1 as reference\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2intersect.pl -f <coordinate file1> -s <coordinate file2>\n";
	print STDERR " -f <coorfile1>   [coordinates in bed format (chr strand start end)]\n";
	print STDERR " -s <coorfile2>   [coordinates in bed format (chr strand start end)]\n\n";
	exit(-1);
}
############################################################################################################

open(INFILE, $coor{'coorfile1'}{'name'}) || die $!;
push(@{$coor{'coorfile1'}{'data'}}, <INFILE>);
close INFILE;
open(INFILE, $coor{'coorfile2'}{'name'}) || die $!;
@{$coor{'coorfile2'}{'data'}}=<INFILE>;
close INFILE;

foreach my $l1(@{$coor{'coorfile1'}{'data'}}) {
	next if($l1=~/^\#/ || $l1=~/^\s*$/);
	($coor{'coorfile1'}{'chr'}, $coor{'coorfile1'}{'strand'}, $coor{'coorfile1'}{'start'}, $coor{'coorfile1'}{'end'}, my @t) = split(/\s+/, $l1);
	if($coor{'coorfile1'}{'chr'}=~/^$/ || $coor{'coorfile1'}{'start'}=~/^$/ || $coor{'coorfile1'}{'end'}=~/^$/) {
		print STDERR "Incorrect coordinate $l1\n";
		usage();
	}

#print "$coor{'coorfile1'}{'chr'}, $coor{'coorfile1'}{'strand'}, $coor{'coorfile1'}{'start'}, $coor{'coorfile1'}{'end'}\n";

	foreach my $l2(@{$coor{'coorfile2'}{'data'}}) {
		next if($l2=~/^\#/ || $l2=~/^\s*$/);
		($coor{'coorfile2'}{'chr'}, $coor{'coorfile2'}{'strand'}, $coor{'coorfile2'}{'start'}, $coor{'coorfile2'}{'end'}, my @t) = split(/\s+/, $l2);
		if($coor{'coorfile2'}{'chr'}=~/^$/ || $coor{'coorfile2'}{'start'}=~/^$/ || $coor{'coorfile2'}{'end'}=~/^$/) {
print "$coor{'coorfile2'}{'chr'}\t$coor{'coorfile2'}{'start'}\t$coor{'coorfile2'}{'end'}\n";
			print STDERR "Incorrect coordinate $l2\n";
			usage();
		}
		next if($coor{'coorfile2'}{'chr'}!~/^$coor{'coorfile1'}{'chr'}$/);
		next if($coor{'coorfile2'}{'strand'}!~/^\Q$coor{'coorfile1'}{'strand'}\E$/);
#print "\t$coor{'coorfile2'}{'chr'}, $coor{'coorfile2'}{'strand'}, $coor{'coorfile2'}{'start'}, $coor{'coorfile2'}{'end'}\n";
	
	if($coor{'coorfile2'}{'start'}>=$coor{'coorfile1'}{'start'} && $coor{'coorfile2'}{'end'}<=$coor{'coorfile1'}{'end'} ||
	$coor{'coorfile2'}{'start'}>=$coor{'coorfile1'}{'start'} && $coor{'coorfile2'}{'start'}<=$coor{'coorfile1'}{'end'} ||
	$coor{'coorfile2'}{'end'}>=$coor{'coorfile1'}{'start'} && $coor{'coorfile2'}{'end'}<=$coor{'coorfile1'}{'end'} ||
	$coor{'coorfile2'}{'start'}<=$coor{'coorfile1'}{'start'} && $coor{'coorfile2'}{'end'}>=$coor{'coorfile1'}{'end'}) {
			chomp($l1); chomp($l2);	
			print "$l1\t$l2\n";
			last;
		}
	}
}

exit;
