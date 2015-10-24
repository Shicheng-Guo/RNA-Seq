#!/usr/bin/perl 

#!/usr/local/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $clusters $annotations $overlap $sizeCutoff $printOpt $flag);


# -----------------------------------------------------------------------------
# OPTIONS

$overlap = -1;
$sizeCutoff = -1;
$printOpt = 0;
$flag=0;

GetOptions ("c=s"       => \$clusters,
	    "a=s"       => \$annotations,
	    "o=s"       => \$overlap,
	    "s=s"       => \$sizeCutoff,
	    "p=s"       => \$printOpt,
	    "f=s"	=> \$flag,
	    "help"      => \$help,
             "h"        => \$help);

usage() if ($help || !$clusters || !$annotations);

if($flag==0) {
	my %ncRNAs = ();

	my $c = 0;
	open(NCRNAS, "<$annotations") || die "cannot open $annotations\n";
	while(<NCRNAS>){
		chomp;
		my @line = split(/\s+/, $_);
		$ncRNAs{$c}{chrom} = $line[0];
		$ncRNAs{$c}{start} = $line[1];
		$ncRNAs{$c}{end} = $line[2];
		$ncRNAs{$c}{id} = $line[3];
		$ncRNAs{$c}{strand} = $line[5];
		$ncRNAs{$c}{source} = $line[6];
		$ncRNAs{$c}{type} = $line[7];
		$ncRNAs{$c}{class} = $line[8];
		$c++;
	}
	close(NCRNAS);

	my $s = 0;
	open(FILE, "<$clusters") || die "cannot open $clusters\n";
	while(<FILE>){
		chomp;
		next if(/\#/);
		if(/^\>/){
			my @line = split(/\s+/, $_);
			if(scalar(@line)!=8) { print STDERR "Incorrect format of Input file\n$_\n"; exit; }
			my $flag = "a"; $s = 0;
			my $perOverlap=0;
			foreach my $c (keys %ncRNAs){
				next if($ncRNAs{$c}{chrom} ne $line[1]);
				next if($ncRNAs{$c}{strand} ne $line[4]);

				if($sizeCutoff != -1){ next if($sizeCutoff * ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1) < ($line[3] - $line[2] + 1)); }
 
				if($line[2] >= $ncRNAs{$c}{start} && $line[2] <= $ncRNAs{$c}{end} ||
				$line[3] >= $ncRNAs{$c}{start} && $line[3] <= $ncRNAs{$c}{end} ||
				$line[2] >= $ncRNAs{$c}{start} && $line[3] <= $ncRNAs{$c}{end} ||
				$line[2] <= $ncRNAs{$c}{start} && $line[3] >= $ncRNAs{$c}{end}) {

					# check fraction of ncRNA overlapping with block group
					my $sum = 0;
					for(my $x = $ncRNAs{$c}{start}; $x <= $ncRNAs{$c}{end}; $x++){
						if($x >= $line[2] && $x <= $line[3]){$sum++;}
					}
					if($overlap != -1){ next if($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1) < $overlap); }
					
					if(($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1)) > $perOverlap) {
						$flag = "$_\t$ncRNAs{$c}{id}\t$ncRNAs{$c}{type}\t$ncRNAs{$c}{class}";
						$perOverlap = ($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1));
						$s = 1;
					}
				}
			}
			# print header for unannotated block group
			if($flag eq "a" && $printOpt == 0){
				print "$_\tn\/a\tn\/a\tn\/a\n";
			}
			# print header for annotated block group
			else {
				print "$flag\n";
			}
		}
		elsif($s == 1) {
			print "$_\n";
		}
		elsif($printOpt == 0) {
			print "$_\n";
		}
	}
	close(FILE);
}
elsif($flag==1) {
	my %loci = ();

	my $c = 0;
	open(LOCI, "<$annotations") || die "cannot open $annotations\n";
	while(<LOCI>){
		chomp;
		my @line = split(/\s+/, $_);
		$loci{$c}{chrom} = $line[0];
		$loci{$c}{start} = $line[1];
		$loci{$c}{end} = $line[2];
		$loci{$c}{id} = $line[3];
		$loci{$c}{strand} = $line[5];
		$loci{$c}{source} = $line[6];
		$loci{$c}{type} = $line[7];
		$loci{$c}{class} = $line[8];
		$c++;
	}
	close(LOCI);

	my $s = 0;
	open(FILE, "<$clusters") || die "cannot open $clusters\n";
	while(<FILE>){
		chomp;
		next if(/\#/);
		if(/^\>/){
			my @line = split(/\s+/, $_);
			if(scalar(@line)!=11) { print STDERR "Input file is not flagged with ncRNA annotation\n"; exit; }
			my $flag = "a"; $s = 0;
			my $perOverlap=0;
			foreach my $c (keys %loci){
				next if($loci{$c}{chrom} ne $line[1]);
				next if($loci{$c}{strand} ne $line[4]);

				if($sizeCutoff != -1){ next if($sizeCutoff * ($loci{$c}{end} - $loci{$c}{start} + 1) < ($line[3] - $line[2] + 1)); }
 
				if($line[2] >= $loci{$c}{start} && $line[2] <= $loci{$c}{end} ||
				$line[3] >= $loci{$c}{start} && $line[3] <= $loci{$c}{end} ||
				$line[2] >= $loci{$c}{start} && $line[3] <= $loci{$c}{end} ||
				$line[2] <= $loci{$c}{start} && $line[3] >= $loci{$c}{end}) {

					# check fraction of loci overlapping with block group
					my $sum = 0;
					for(my $x = $loci{$c}{start}; $x <= $loci{$c}{end}; $x++){
						if($x >= $line[2] && $x <= $line[3]){$sum++;}
					}
					if($overlap != -1){ next if($sum / ($loci{$c}{end} - $loci{$c}{start} + 1) < $overlap); }
					
					if(($line[10]=~/\./ || $line[10]=~/n\/a/) && ($sum / ($loci{$c}{end} - $loci{$c}{start} + 1)) > $perOverlap) {
						$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$loci{$c}{'type'}";
						$s = 1;
					}
					elsif(($loci{$c}{end} - $loci{$c}{start} + 1) > $perOverlap) {
						$line[10] = $line[10]."_".$loci{$c}{'type'};
						$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]";
						$s = 1;
					}
				}
			}
			# print header for unannotated block group
			if($flag eq "a" && $printOpt == 0){
				print "$_\n";
			}
			# print header for annotated block group
			else {
				print "$flag\n";
			}
		}
		elsif($s == 1) {
			print "$_\n";
		}
		elsif($printOpt == 0) {
			print "$_\n";
		}
	}
	close(FILE);
}

sub usage {
	print STDERR "\nusage: flagKnownClusters.pl -c <file> -a <file> [OPTIONS]\n";
	print STDERR "flag clusters overlapping with annotated ncRNAs or loci\n";
	print STDERR "\n";
	print STDERR "[INPUT]\n";
	print STDERR " -c <file>    cluster file\n";
	print STDERR " -a <file>    ncRNAs or loci annotation file\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <file>    min overlap in % (default = off)\n";
	print STDERR " -s <file>    size cutoff, i.e. 1.5 times of the length (default = off)\n";
	print STDERR " -p <int>     print option (0 = all [default]; 1 = only clusters with annotation)\n";
	print STDERR " -f <int>     flag (0 = annotation [default]; 1 = loci)\n";
	print STDERR "[VERSION]\n";
	print STDERR " 11-02-2011\n";
	print STDERR "[BUGS]\n";
	print STDERR " Please report bugs to sachin\@rth.dk\n";
	print STDERR "\n";
	exit(-1);
}
