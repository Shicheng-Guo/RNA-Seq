#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

##################################################################################################
use vars qw($mafFile $refGenome $onlyMaxLen $stitch $computeSS $help);
my %threshold=();
$refGenome="hg19";

GetOptions ("i=s"  => \$mafFile,
            "r=s"  => \$refGenome,
            "l=i"  => \$threshold{'len'},
            "m"    => \$onlyMaxLen,
            "d=i"  => \$threshold{'dist'},
            "s"    => \$computeSS,
            "t"    => \$stitch,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

##################################################################################################
sub usage {
	print STDERR "\nProgram: editMaf.pl (compute average pairwise identify for MAF blocks)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: editMaf.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>        [maf file | STDIN]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -r <string>      [reference genome (default: hg19)]\n";
	print STDERR " -l <int>         [minimum length of maf blocks (default: 1)]\n";
	print STDERR " -m               [compute only for maximum length MAF block]\n";
	print STDERR " -s               [compute secondary structure]\n";
	#print STDERR " -d <int>         [maximum distance between stitched blocks (default: 10)]\n";
	#print STDERR " -t               [stitch maf blocks]\n";
	print STDERR " -h <help>\n";
	print STDERR "[TIP]\n";
	print STDERR " before using PETfold for RNA secondary structure prediction, please export following path\n";
	print STDERR " export PETFOLDBIN=/home/users/seemann/projects/PETfold/code/c/bin/\n\n";
	exit(-1);
}

##################################################################################################

## open coordinate file in BED format
my @data=();
if(defined($mafFile)) {
	@data=openFile($mafFile);
}
else {
	@data=openFile();
}

my %genome=();
## read MAF blocks
my %mafInfo=(); my $MAF_COUNTER=0; my @comment=();
foreach my $l(@data) {
	chomp($l);
	if($l=~/^\#/) { push(@comment, $l); }
	elsif($l=~/^a score/) {
		$MAF_COUNTER++;
		@{$mafInfo{$MAF_COUNTER}{'comment'}}=@comment;
		@comment=();
		$mafInfo{$MAF_COUNTER}{'score'}=$l;
	}
	elsif($l=~/^s\s+/) {
		my @t=split(/\s+/, $l);
		my($genome, $chr) = split(/\./, $t[1]);
		$mafInfo{$MAF_COUNTER}{'chr'}{$genome}=$chr;
		$mafInfo{$MAF_COUNTER}{'start'}{$genome}=$t[2];
		$mafInfo{$MAF_COUNTER}{'end'}{$genome}=$t[2]+($t[3]-1);
		$mafInfo{$MAF_COUNTER}{'len'}{$genome}=$t[3];
		$mafInfo{$MAF_COUNTER}{'strand'}{$genome}=$t[4];
		$mafInfo{$MAF_COUNTER}{'srcSize'}{$genome}=$t[5];
		$mafInfo{$MAF_COUNTER}{'aln'}{$genome}="$t[6]";
		push(@{$genome{'all'}}, $genome);
	}
}

if($stitch) {
	## determine common organisms shared between all MAF blocks
	@{$genome{'common'}}=commonHashKeys(\@{$genome{'all'}}, scalar(keys(%mafInfo)));

	#foreach(@{$genome{'common'}}) { print "$_\n"; }

	## stitch MAF blocks
	foreach my $gen(@{$genome{'common'}}) {
		my $i=0; my $last_end=(); my %stitchMafInfo=();
		foreach my $maf(sort { $mafInfo{$a}{'start'}{$refGenome} <=> $mafInfo{$b}{'start'}{$refGenome} } keys(%mafInfo)) {
			if($i==0) {
				$stitchMafInfo{$i}{'chr'}=$mafInfo{$maf}{'chr'}{$gen};
				$stitchMafInfo{$i}{'len'}+=$mafInfo{$maf}{'len'}{$gen};
				$stitchMafInfo{$i}{'strand'}=$mafInfo{$maf}{'strand'}{$gen};
				$stitchMafInfo{$i}{'srcSize'}=$mafInfo{$maf}{'srcSize'}{$gen};
				$stitchMafInfo{$i}{'start'}{$gen}=$mafInfo{$maf}{'start'}{$gen};
				$stitchMafInfo{$i}{'aln'}{$gen}.=$mafInfo{$maf}{'aln'}{$gen};
				$last_end=$mafInfo{$maf}{'end'}{$gen};
			}
			else {
				for(my $j=1; $j<($mafInfo{$maf}{'start'}{$gen}-$last_end); $j++) {
					$stitchMafInfo{$i}{'aln'}{$gen}.="-";
				}
				$stitchMafInfo{$i}{'aln'}{$gen}.=$mafInfo{$maf}{'aln'}{$gen};
				$last_end=$mafInfo{$maf}{'end'}{$gen};
			}
			$i++;
		}
	}
}
else {
	## determine maximum MAF length
	my $MAX_LENGTH=0;
	foreach my $maf(keys(%mafInfo)) {
		if($mafInfo{$maf}{'len'}{$refGenome}>$MAX_LENGTH) {
			$MAX_LENGTH=$mafInfo{$maf}{'len'}{$refGenome};
		}
	}

	## set length threshold to maximum length, if $threshold{'len'} is not defined
	$threshold{'len'}=$MAX_LENGTH if(!defined($threshold{'len'}));

	my @avgPercentIdentity=(); my @avgMAFLength=(); my @avgGenomeCount=(); my @SS=();
	foreach my $maf(sort { $mafInfo{$a}{'start'}{$refGenome} <=> $mafInfo{$b}{'start'}{$refGenome} } keys(%mafInfo)) {
		if($mafInfo{$maf}{'len'}{$refGenome}>=$threshold{'len'}) {
			## open temporary file to keep fasta sequence per MAF block
			my $tmpFileName=$refGenome."_".$mafInfo{$maf}{'chr'}{$refGenome}."_".$mafInfo{$maf}{'start'}{$refGenome}."_".$mafInfo{$maf}{'len'}{$refGenome};
			open(TMPFILE, ">$tmpFileName.fasta") || die $!;

			## MAF block length
			push(@avgMAFLength, $mafInfo{$maf}{'len'}{$refGenome});

			foreach my $gen(keys(%{$mafInfo{$maf}{'aln'}})) {
				print TMPFILE ">$gen"."_"."$mafInfo{$maf}{'chr'}{$gen}"."_"."$mafInfo{$maf}{'start'}{$gen}"."_"."$mafInfo{$maf}{'len'}{$gen}\n";
				print TMPFILE "$mafInfo{$maf}{'aln'}{$gen}\n";
			}
			close(TMPFILE);

			## compute average percent sequence identity in MAF block
			push(@avgPercentIdentity, `seqIdentity.pl -f $tmpFileName.fasta | tail -n 1`);

			## number of genomes included in the MAF block
			push(@avgGenomeCount, scalar(keys(%{$mafInfo{$maf}{'aln'}})));

			## compute secondary structure
			if($computeSS && $MAX_LENGTH>=50) {
				#system("export PETFOLDBIN=/home/users/seemann/projects/PETfold/code/c/bin/");
				#push(@SS, `PETfold -f $tmpFileName.fasta -g 0.50`);
				system("fasta2clustal.pl -i $tmpFileName.fasta -o $tmpFileName.clustalw -f clustalw > /dev/null");
				push(@SS, `RNAalifold -p -r -d2 -noLP -color -aln < $tmpFileName.clustalw`);

				## remove temporary file
				system("rm $tmpFileName.clustalw");
				system("rm alifold.out alidot.ps aln.ps");
				system("rm $tmpFileName.clustalw");
				system("mv alirna.ps $tmpFileName.ps");
			}

			## remove temporary file
			system("rm $tmpFileName.fasta");
		}
	}

	## print results
	if(defined($mafFile)) {	print "$mafFile\t"; }
	else { print "STDIN\t"; }

	## compute average of @avgPercentIdentity
	my $sum=0;
	foreach(@avgPercentIdentity) { $sum+=$_; }
	if($sum>0) {
		printf("%0.2f\t", $sum/scalar(@avgPercentIdentity));
	}
	else { print "0\t"; }

	## compute average of @avgMAFLength
	$sum=0;
	foreach(@avgMAFLength) { $sum+=$_; }
	if($sum>0) {
		printf("%0.2f\t", $sum/scalar(@avgMAFLength));
	}
	else { print "0\t"; }

	## compute average of @avgGenomeCount
	$sum=0;
	foreach(@avgGenomeCount) { $sum+=$_; }
	if($sum>0) {
		printf("%0.2f\n", $sum/scalar(@avgGenomeCount));
	}
	else { print "0\n"; }

	## print secondary structure information
	if($computeSS && $MAX_LENGTH>=50) {
		foreach(@SS) { print "#$_";	}
	}
}
exit(0);
