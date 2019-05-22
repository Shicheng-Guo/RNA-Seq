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

use vars qw ($help $bgFile $annoFile $overlap $sizeCutoff $printOpt $flag);

##########################################################################################################
## parse input options
$overlap = -1;
$sizeCutoff = -1;
$printOpt = 0;
$flag=0;

GetOptions ("c=s"  => \$bgFile,
            "a=s"  => \$annoFile,
            "o=s"  => \$overlap,
            "s=s"  => \$sizeCutoff,
            "p=s"  => \$printOpt,
            "f=s"	 => \$flag,
            "help" => \$help,
            "h"    => \$help);

usage() if ($help || !$bgFile || !$annoFile);

##########################################################################################################

sub usage {
	print STDERR "\nusage: flagCluster.pl -c <file> -a <file> [OPTIONS]\n";
	print STDERR "flag block groups overlapping with annotated ncRNAs or loci\n";
	print STDERR "\n";
	print STDERR "[INPUT]\n";
	print STDERR " -c <file>    [block group file in blockbuster format]\n";
	print STDERR " -a <file>    [ncRNAs or loci annotation file(s) seperated by comma]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <file>    [min overlap in % (default: off)]\n";
	print STDERR " -s <file>    [size cutoff, i.e. 1.5 times of the length (default: off)]\n";
	print STDERR " -p <int>     [print option, 0 = all, 1 = only block groups with annotation (default: 0)]\n";
	print STDERR " -f <int>     [flag, 0 = annotation, 1 = loci (default: 0)]\n";
	print STDERR "[VERSION]\n";
	print STDERR " 18-04-2012\n";
	print STDERR "[BUGS]\n";
	print STDERR " Please report bugs to sachin\@rth.dk\n";
	print STDERR "\n";
	exit(-1);
}

##########################################################################################################

if($flag==0) {
	my %ncRNAs = ();

	my $c = 0;
	foreach my $annotations(split(/\,/, $annoFile)) {
		if($annotations=~/\.gz$/) { open(NCRNAS, "gunzip -c $annotations |") || die $!; }
		else { open(NCRNAS, "<$annotations") || die "cannot open $annotations\n"; }
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
	}

	my $s = 0;
	if($bgFile=~/\.gz$/) { open(FILE, "gunzip -c $bgFile |") || die $!; }
	else { open(FILE, "<$bgFile") || die "cannot open $bgFile\n"; }
	while(<FILE>){
		chomp;
		next if(/\#/);
		my @line = split(/\s+/, $_); my $doFlag=0;
		if(/^\>/ && scalar(@line)==8) { $doFlag=1; }
		elsif(/^\>/ && scalar(@line)==11 && $line[9]=~/n\/a/) { $doFlag=1; }
		elsif(/^\>/ && scalar(@line)!=11 && scalar(@line)!=8) { print STDERR "Incorrect input format of $bgFile"; exit; }

		if(/^\>/ && $doFlag){
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
					if($overlap != -1){ next if((($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1))*100) < $overlap); }
					
					if((($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1))*100) > $perOverlap) {
						$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$ncRNAs{$c}{id}\t$ncRNAs{$c}{type}\t$ncRNAs{$c}{class}";
						$perOverlap = ($sum / ($ncRNAs{$c}{end} - $ncRNAs{$c}{start} + 1));
						$s = 1;
					}
				}
			}
			# print header for unannotated block group
			if($flag eq "a" && $printOpt == 0){
				print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tn\/a\tn\/a\t.\n";
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
	my %loci = (); my %anno_chrom=();

	my $c = 0;
	foreach my $annotations(split(/\,/, $annoFile)) {
		if($annotations=~/\.gz$/) { open(LOCI, "gunzip -c $annotations |") || die $!; }
		else { open(LOCI, "<$annotations") || die "cannot open $annotations\n"; }
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

			if(!defined($anno_chrom{$loci{$c}{chrom}})) {
				$anno_chrom{$loci{$c}{chrom}}{'start'} = $c;
				$anno_chrom{$loci{$c}{chrom}}{'end'} = $c;
			}
			else {
				$anno_chrom{$loci{$c}{chrom}}{'end'} = $c;
			}
			$c++;
		}
		close(LOCI);
	}

	#foreach(keys(%anno_chrom)) { print "$_\t$anno_chrom{$_}{'start'}\t$anno_chrom{$_}{'end'}\n"; } exit;

	my %perOverlap=();
	my %seen=(); foreach(keys(%loci)) { $seen{$loci{$_}{type}}++; }
	foreach(keys(%seen)) { $perOverlap{$_}=0; } %seen=();

	my $s = 0;
	if($bgFile=~/\.gz$/) { open(FILE, "gunzip -c $bgFile |") || die $!; }
	else { open(FILE, "<$bgFile") || die "cannot open $bgFile\n"; }
	while(<FILE>){
		chomp;
		next if(/\#/);
		if(/^\>/){
			my @line = split(/\s+/, $_);
			if(scalar(@line)!=11) { print STDERR "Input file is not flagged with ncRNA annotation\nLine: $_\n"; exit; }
			my $flag = "a"; $s = 0;
			foreach my $key(keys(%perOverlap)) { $perOverlap{$key}=0; }
			
			my $start=$anno_chrom{$line[1]}{'start'}; my $end=$anno_chrom{$line[1]}{'end'};
			while($end-$start > 2000) {
				my $absPos = sprintf "%.0f", $start+(($end-$start)/2); 
				#print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$start\t$end\t$absPos\n";
				if($loci{$absPos}{'start'} > $line[3]) { $end = $absPos; }
				elsif($loci{$absPos}{'end'} < $line[2]) { $start = $absPos; }
				else { $start = $absPos-1000; $end = $absPos+1000; }
			}
			#print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$loci{$start}{start}\t$loci{$end}{start}\n";

			#foreach my $c (keys %loci){
			#for(my $c=$anno_chrom{$line[1]}{'start'}; $c<=$anno_chrom{$line[1]}{'end'}; $c++) {
			for(my $c=$start; $c<=$end; $c++) {
				next if($loci{$c}{chrom} ne $line[1]);
				next if($loci{$c}{strand} ne $line[4]);
				next if($c < $anno_chrom{$line[1]}{'start'} || $c > $anno_chrom{$line[1]}{'end'});

				if($sizeCutoff != -1){ next if($sizeCutoff * ($loci{$c}{end} - $loci{$c}{start} + 1) < ($line[3] - $line[2] + 1)); }
 
				if($line[2] >= $loci{$c}{start} && $line[2] <= $loci{$c}{end} ||
				$line[3] >= $loci{$c}{start} && $line[3] <= $loci{$c}{end} ||
				$line[2] >= $loci{$c}{start} && $line[3] <= $loci{$c}{end} ||
				$line[2] <= $loci{$c}{start} && $line[3] >= $loci{$c}{end}) {

					# check fraction of loci overlapping with block group
					my $sum = 0;
					for(my $x = $line[2]; $x <= $line[3]; $x++) {
						if($x >= $loci{$c}{start} && $x <= $loci{$c}{end}){$sum++;}
					}
					$perOverlap{$loci{$c}{type}} += ($sum / ($line[3] - $line[2] + 1))*100;
				}
			}
			my @maxOverlapLoci = reverse sort { $perOverlap{$a} <=> $perOverlap{$b} } keys(%perOverlap);
			#print "$maxOverlapLoci[0]\t$perOverlap{$maxOverlapLoci[0]}\n";

			if(($line[10]=~/\./ || $line[10]=~/n\/a/) && $perOverlap{$maxOverlapLoci[0]} > $overlap) {
				$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$maxOverlapLoci[0]";
				$s = 1;
			}
			elsif($perOverlap{$maxOverlapLoci[0]} > $overlap) {
				$line[10] = $maxOverlapLoci[0]."_".$line[10];
				$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]";
				$s = 1;
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

