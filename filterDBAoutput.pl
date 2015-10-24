#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#########################################################################################
# parse input options
use vars qw($dbaFile $bgId $threshold $reportAlignments $onlyHits $help);

# threshold for dba hit
$threshold = 0;
# number of alignments to report
$reportAlignments = "all";

GetOptions ("i=s" => \$dbaFile,
            "b=s"  => \$bgId,
            "t=s"  => \$threshold,
            "n=s"  => \$reportAlignments,
            "x"    => \$onlyHits,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

#########################################################################################
sub usage {
	print STDERR "\nProgram: filterDBAoutput.pl (filter deepBlockAlign output for significant alignments)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: filterDBAoutput.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>   [deepBlockAlign output file | STDIN]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -b <string> [block group id, if filter for a specific block group]\n";
	print STDERR " -t <float>  [minimum score to filter the aligment as significant (default: 0)]\n";
	print STDERR " -n <int>    [number of alignments to report for each block group (default: all)]\n";
	print STDERR " -x          [report only significant alignments]\n";
	print STDERR " -h <help>   [print this helpful message]\n\n";
  exit(-1);
}

#########################################################################################

## open DBA output file
my $INFILE=();
if(defined($dbaFile)) {
	if($dbaFile=~/\.gz$/) { open($INFILE, "gunzip -c $dbaFile |") || die $!; }
	else { open($INFILE, $dbaFile) || die $!; }
}
else { $INFILE = *STDIN; }
my @data = <$INFILE>;
close $INFILE;

## read block group alignments
my %alnInfo=(); my $qFile=(); my $sFile=();
foreach my $l(@data) {
	my @t=split(/\s+/, $l);
	if($l=~/^\#\s+query\:/) {
		$qFile=$t[2];
	}
	elsif($l=~/^\#\s+subject\:/) {
		$sFile=$t[2];
	}
	elsif($l=~/^\#/ || $l=~/^$/) { next; }
	elsif(defined($bgId) && $t[0]=~/\>$bgId\|/) {
		$alnInfo{$t[0]}{$t[1]}{'score'}=sprintf("%0.02f",$t[2]);
		$alnInfo{$t[0]}{$t[1]}{'qFile'}=$qFile;
		$alnInfo{$t[0]}{$t[1]}{'sFile'}=$sFile;
	}
	elsif(!defined($bgId)) {
		$alnInfo{$t[0]}{$t[1]}{'score'}=sprintf("%0.2f",$t[2]);
		$alnInfo{$t[0]}{$t[1]}{'qFile'}=$qFile;
		$alnInfo{$t[0]}{$t[1]}{'sFile'}=$sFile;
	}
}

## print filtered alignments
foreach my $qBG(keys(%alnInfo)) {
	my $count=1;
	foreach my $sBG(reverse sort { $alnInfo{$qBG}{$a}{'score'} <=> $alnInfo{$qBG}{$b}{'score'} } keys(%{$alnInfo{$qBG}})) {
		if($reportAlignments!~/all/ && $alnInfo{$qBG}{$sBG}{'score'}>=$threshold && $count<=$reportAlignments) {
			print "$qBG\t$sBG\t$alnInfo{$qBG}{$sBG}{'score'}\t$alnInfo{$qBG}{$sBG}{'qFile'}\t$alnInfo{$qBG}{$sBG}{'sFile'}\n";
			$count++;
		}
		elsif($reportAlignments=~/all/ && $alnInfo{$qBG}{$sBG}{'score'}>=$threshold) {
			print "$qBG\t$sBG\t$alnInfo{$qBG}{$sBG}{'score'}\t$alnInfo{$qBG}{$sBG}{'qFile'}\t$alnInfo{$qBG}{$sBG}{'sFile'}\n";
			$count++;
		}
		else { last; }
	}
	## not hit observed for $qBG
	if($count==1 && !defined($onlyHits)) {
		print "$qBG\t-\t-\t-\t-\n";
	}
}
exit;
