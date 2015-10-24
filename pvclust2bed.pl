#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

######################################################################################################################
#Parse input options
use vars qw($pvclustfile $stdin $clusterId $outToFile $help);
$clusterId=".*";

GetOptions ("p=s"  => \$pvclustfile,
            "i=s"  => \$clusterId,
						"o"    => \$outToFile,
      	    "help" => \$help,
	          "h"    => \$help);

usage() if($help);

######################################################################################################################
sub usage {
	print STDERR "\nProgram: pvclust2bed.pl (convert pvclust format to BED format)\n";
	print STDERR "Author: University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: pvclust2bed.pl -p <file|STDIN> [OPTIONS]\n";
	print STDERR " -p <file>        [pvclust file|STDIN]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -i <string>      [cluster id, if convert to BED for a specific cluster (default: all)]\n";
	print STDERR " -o               [output to seperate file for each cluster (default: STDOUT)]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}
######################################################################################################################

## open pvclust file
my @data=();
if(defined($pvclustfile)) {
	@data=openFile($pvclustfile);
	$pvclustfile=~s/^.+\///g;
}
else {
	@data=openFile();
	$pvclustfile="pvclust";
}

my $OUTFILE=(); my $CLUSTER_ID=(); 
foreach(@data) {
	chomp($_);
	next if($_=~/\#\s+pvclust/);
	if($_=~/^\$clusters\[\[[0-9]+\]\]$/) {
		close $OUTFILE if(defined($OUTFILE));
		$CLUSTER_ID=();

		$_=~s/[\$\[\]]+//g;
		if($_=~/$clusterId/ && defined($outToFile)) {
			my $bedfile=$pvclustfile.".".$_.".bed";
			open($OUTFILE, ">$bedfile") || die $!;
			print "$bedfile\n";
			$CLUSTER_ID=$_;
		}
		elsif($_=~/$clusterId/) {
			$CLUSTER_ID=$_;
		}
	}
	elsif($_=~/\s*\[[0-9]+\]\s+\"/ && defined($CLUSTER_ID)) {
		my @t=split(/\./, $_);
		my $chr = $t[4];
		my ($start, $end) = split(/\-/, $t[5]);
		$t[0]=~s/^.*\[[0-9]+\].*\>//g;
		if(defined($OUTFILE)) {
			print $OUTFILE "$chr\t$start\t$end\t$CLUSTER_ID\t1\t$t[6]\n";
		}
		else {
			print "$chr\t$start\t$end\t$CLUSTER_ID\t1\t$t[6]\n";
		}
	}
}

close $OUTFILE if(defined($OUTFILE));
exit(0);
