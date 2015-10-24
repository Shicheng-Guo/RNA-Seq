#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Tie::File;
use perlModule qw(checkMapFormat);

##################################################################################################
# parse input options
use vars qw($help $mirbaseFile $mirCoorFile);

GetOptions ("i=s"  => \$mirbaseFile,
            "j=s"  => \$mirCoorFile,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$mirbaseFile || !$mirCoorFile);

#################################################################################################
sub usage {
	print STDERR "\nProgram: mirbase2bed.pl (convert deep sequencing data from miRBase into BED format)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: mirbase2bed.pl -i <file> -j <file> [OPTIONS]\n";
	print STDERR " -i <file>    [miRBase's mirna_read_count_by_experiment.txt.gz file]\n";
	print STDERR " -j <file>    [genomic coordinates of miRNA from miRBase's gff files]\n";
	print STDERR "\n";
	exit(-1);
}
##################################################################################################

## retrieve all unique miRNA accession
my @accession=`zless $mirbaseFile | cut -f 1 | sort | uniq`;
#print scalar(@{$miRNA{'accession'}})."\n";

## retrieve genomic coordinates of each miRNA accession
my %miRNA=();
foreach my $l(@accession) {
	chomp($l);
	my @t=split(/\s+/, `zless $mirCoorFile | grep "_$l"`);
	if(scalar(@t)>0) {
		$miRNA{$l}{'chr'}=$t[0];
		$miRNA{$l}{'start'}=$t[1];
		$miRNA{$l}{'end'}=$t[2];
		$miRNA{$l}{'strand'}=$t[5];
		($miRNA{$l}{'name'}, $miRNA{$l}{'accession'})=split(/\_/,$t[3]);
		if(!defined($miRNA{$l}{'name'}) || !defined($miRNA{$l}{'accession'})) {
			print "No annotation found for $l\n";
		}
	}
	#($miRNA{$l}{'chr'}, $miRNA{$l}{'start'}, $miRNA{$l}{'end'}, $miRNA{$l}{'strand'}, $miRNA{$l}{'accession'}, $miRNA{$l}{'name'})=split(/\s+/, `zless $mirCoorFile | grep -w $l`);
	#print "$l\t$miRNA{$l}{'chr'}\t$miRNA{$l}{'start'}\t$miRNA{$l}{'end'}\t$miRNA{$l}{'strand'}\t$miRNA{$l}{'accession'}\t$miRNA{$l}{'name'}\n";
}

#print keys(%miRNA)."\n";

## convert miRNA to BED format
my @data=`zless $mirbaseFile`;

my %read=();
foreach my $l(@data) {
	my @t=split(/\s+/, $l);
	## create read label
	if(defined($read{$t[3]})) {
		$read{$t[3]}{'counter'}++;
		$read{$t[3]}{'name'}="$t[3]_r".($read{$t[3]}{'counter'})."_1";
	}
	else {
		$read{$t[3]}{'counter'}=1;
		$read{$t[3]}{'name'}="$t[3]_r".($read{$t[3]}{'counter'})."_1";
	}
	## read length
	$read{$t[3]}{'len'}=scalar(@{[split(//,$t[2])]});

	## mapping information, if coordinates of the miRNA are found in coordinate file
	if(defined($miRNA{$t[0]}{'name'})) {
		$read{$t[3]}{'name'}=$miRNA{$t[0]}{'name'}."_"."$read{$t[3]}{'name'}";
		$read{$t[3]}{'start'}=($miRNA{$t[0]}{'start'}-1)+$t[7];
		$read{$t[3]}{'end'}=($miRNA{$t[0]}{'start'}-1)+($t[7]-1)+$read{$t[3]}{'len'};

		#print "$l\t$miRNA{$t[0]}{'chr'}\t$miRNA{$t[0]}{'start'}\t$miRNA{$t[0]}{'end'}\t$miRNA{$t[0]}{'strand'}\t$miRNA{$t[0]}{'accession'}\t$miRNA{$t[0]}{'name'}\n";
		print "$miRNA{$t[0]}{'chr'}\t$read{$t[3]}{'start'}\t$read{$t[3]}{'end'}\t$read{$t[3]}{'name'}\t$t[4]\t$miRNA{$t[0]}{'strand'}\n";
	}
	else { print STDERR "No coordinate found for $t[0]\n"; }
}

exit;

