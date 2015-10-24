#!/usr/bin/perl -w

use strict;
use warnings;
use perlModule;
use Getopt::Long;

############################################################################################################
## parse input options

use vars qw($help);

GetOptions ("help" => \$help,
            "h"    => \$help);

usage() if($help);

############################################################################################################
sub usage {
	print STDERR "\nProgram: uniqLongestGene.pl (given the UCSC genes in BED format, retrieve only unique and longest transcript of each gene)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: uniqLongestGene.pl\n";
	print STDERR " -h <help>\n";
	exit(-1);
}
############################################################################################################
my %gene=(); my $key=(); my $prevKey=();
while(<STDIN>) {
	my @F=split(/\t{1}/, $_);
	$key="$F[0]_$F[3]_$F[5]"; 
	chomp($F[7]); chomp($F[8]);
	my $overlap=checkOverlap($F[1], $F[2], $gene{$key}{'start'}, $gene{$key}{'end'}, 0);
	if(defined($gene{$key}) && $overlap>=1) {
		if($F[1]<$gene{$key}{'start'}) { $gene{$key}{'start'}=$F[1]; }
		if($F[2]>$gene{$key}{'end'}) { $gene{$key}{'end'}=$F[2]; }
		if($F[6]<$gene{$key}{'cdsStart'}) { $gene{$key}{'cdsStart'}=$F[6]; }
		if($F[7]>$gene{$key}{'cdsEnd'}) { $gene{$key}{'cdsEnd'}=$F[7]; }
		if($F[8]!~/^$/) { $gene{$key}{'protein'}=$F[8]; }
		$prevKey=$key;
		#print "defined: $F[1]_$F[2]_$key\n";
	}
	else {
		#print "new: $F[1]_$F[2]_$key\n";
		if(defined($prevKey)) {
			print "$gene{$prevKey}{'chr'}\t$gene{$prevKey}{'start'}\t$gene{$prevKey}{'end'}\t$gene{$prevKey}{'name'}\t0\t$gene{$prevKey}{'strand'}\t$gene{$prevKey}{'cdsStart'}\t$gene{$prevKey}{'cdsEnd'}\t$gene{$prevKey}{'protein'}\t$prevKey\n";
		}
		%gene=();
		$gene{$key}{'chr'}=$F[0];
		$gene{$key}{'start'}=$F[1];
		$gene{$key}{'end'}=$F[2];
		$gene{$key}{'strand'}=$F[5];
		$gene{$key}{'name'}=$F[3];
		$gene{$key}{'cdsStart'}=$F[6];
		$gene{$key}{'cdsEnd'}=$F[7];
		$gene{$key}{'protein'}=$F[8];
		$prevKey=$key;
	}
}
if(defined($prevKey)) {
	print "$gene{$prevKey}{'chr'}\t$gene{$prevKey}{'start'}\t$gene{$prevKey}{'end'}\t$gene{$prevKey}{'name'}\t0\t$gene{$prevKey}{'strand'}\t$gene{$prevKey}{'cdsStart'}\t$gene{$prevKey}{'cdsEnd'}\t$gene{$prevKey}{'protein'}\t$prevKey\n"; $prevKey=$key;
}

exit(0);
