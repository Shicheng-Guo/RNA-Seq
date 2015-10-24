#!/usr/bin/perl -w

use strict;
use warnings;

use GD::Simple;
use Getopt::Long;
use perlModule;
use POSIX qw(strftime);

#########################################################################################################
## parse input options
use vars qw($coor $trackFile $strand $expr $help);

GetOptions ("c=s"  => \$coor,
	          "t=s"  => \$trackFile,
						"s=s"  => \$strand,
            "e"    => \$expr,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$coor || !$trackFile);

#########################################################################################################
sub usage {
	print STDERR "\nProgram: coor2wig.pl (coordinate to read map density in wig format)\n";
	print STDERR "Author: University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2wig.pl -c <string> -t <file> [OPTIONS]\n";
	print STDERR " -c <string>   [coordinate as chr:start-end]\n";
	print STDERR " -t <file>     [UCSC track file (chr start end expr)]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -s <string>   [strand]\n";
	print STDERR " -e            [measure read count (default: expression per nucleotide)]\n";
	print STDERR " -h            [help message]\n";
	exit(-1);
}
#########################################################################################################

## read track file
my @data=openFile($trackFile);

## read input coordinate
my @coor=split(/[\:\-]+/,$coor);

## store track data corresponding to input coordinate
my @track=(); my $found=();
foreach my $l(@data) {
	next if($l=~/^\#/);
	my @t=split(/\s+/,$l);
	if($l=~/track type\=wiggle/) {
		$found=1;
	}
	elsif($l=~/track/) {
		$found=0;
	}
	elsif($found && $t[0]=~/^$coor[0]$/) {
		push(@track, $l);
	}
}

if($expr) {
	my $expr=0;
	my %prevCoor=();
	foreach my $l(@track) {
		chomp($l);
		my @t=split(/\s+/,$l);
		next if($t[0]!~/^$coor[0]$/);
		my $overlap=checkOverlap($coor[1], $coor[2], $t[1], $t[2], 0);
		if($overlap>=1) {
print "$l\n";
			if(defined($prevCoor{'start'}) && defined($prevCoor{'end'})) {
				my $overlap=checkOverlap($prevCoor{'start'}, $prevCoor{'end'}, $t[1], $t[2], 0);
				if($overlap>=1 && $t[3]>$prevCoor{'expr'}) {
					$expr+=($t[3]-$prevCoor{'expr'});
				}
				elsif($overlap<1) {
					$expr+=$t[3];
				}
			}
			else {
				$expr+=$t[3];
			}
			#print "$l\t$expr\n";
			$prevCoor{'start'}=$t[1];
			$prevCoor{'end'}=$t[2];
			$prevCoor{'expr'}=$t[3];
		}
	}
	print "$trackFile\t$coor\t$expr\n";
}
else {
	my $expr=0;
	foreach my $l(@track) {
		chomp($l);
		my @t=split(/\s+/,$l);
		next if($t[0]!~/^$coor[0]$/);
		my $overlap=checkOverlap($coor[1], $coor[2], $t[1], $t[2], 0);
		if($overlap>=1) {
			$expr+=($t[2]-$t[1])*$t[3];
		}
	}
	printf("%s\t%s\t%0.2f\n", $trackFile, $coor, $expr/(($coor[2]-$coor[1])+1));
}

exit;
=cu
	## compute expression at each nucleotide position
	for(my $i=$coor[1]; $i<=$coor[2]; $i++) {
		my $expr=0;
		foreach my $l(@track) {
			my @t=split(/\s+/,$l);
			next if($t[0]!~/^$coor[0]$/);
			if($i>=$t[1] && $i<$t[2]) {
				$expr+=$t[3];
			}
		}
		#print "$i\t$expr\n" if($expr>0);
		print "$i\t$expr\n";
	}
=cut
