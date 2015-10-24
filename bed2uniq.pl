#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

#############################################################################################################
use vars qw($bedFile $overlapThreshold $help);
$overlapThreshold=1;

GetOptions ("f=s"  => \$bedFile,
            "o=i"  => \$overlapThreshold,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: bed2uniq.pl (identify unique coordinates from a BED file)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: bed2uniq.pl -f <file | STDIN> [OPTIONS]\n";
	print STDERR " -f <file>        [BED file (optional)]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <overlap>     [default 1 nt (optional)]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}

#############################################################################################################

## open bed file
my %bed=();
if(defined($bedFile)) {
	@{$bed{'all'}}=openFile($bedFile);
}
else {
	@{$bed{'all'}}=openFile();
}

## retrieve unique coordinates
my %coor=();
for(my $i=0; $i<scalar(@{$bed{'all'}}); $i++) {
	chomp($bed{'all'}[$i]);
	my @t=split(/\s+/,$bed{'all'}[$i]);
	$coor{'first'}{'chr'}=$t[0];
	$coor{'first'}{'start'}=$t[1];
	$coor{'first'}{'end'}=$t[2];
	push(@{$bed{'uniq'}}, $bed{'all'}[$i]);
	for(my $j=0; $j<scalar(@{$bed{'all'}}); $j++) {
		chomp($bed{'all'}[$j]);
		next if($i==$j);
		my @t=split(/\s+/,$bed{'all'}[$j]);
		$coor{'second'}{'chr'}=$t[0];
		$coor{'second'}{'start'}=$t[1];
		$coor{'second'}{'end'}=$t[2];
		my $overlap=checkOverlap($coor{'first'}{'start'}, $coor{'first'}{'end'}, $coor{'second'}{'start'}, $coor{'second'}{'end'}, 0, $coor{'first'}{'chr'}, $coor{'second'}{'chr'});
		if($overlap>=1) {
			splice(@{$bed{'all'}}, $j, 1); 
		}
	}
}

## print unique coordinates
foreach(@{$bed{'uniq'}}) {
	print "$_\n";
}

exit(0);
