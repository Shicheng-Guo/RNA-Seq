#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;
use Tie::File;

##################################################################################################
use vars qw($mafFile $refGenome $computeSS $computeAvg $help);
my %threshold=();
$refGenome="hg19";

GetOptions ("i=s"  => \$mafFile,
            "r=s"  => \$refGenome,
            "l=i"  => \$threshold{'len'},
            "d=i"  => \$threshold{'dist'},
            "s=i"  => \$computeSS,
            "a"    => \$computeAvg,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

##################################################################################################
sub usage {
	print STDERR "\nProgram: editMaf.pl (compute MAF block properties (average pairwise identity, length, genomes and SS)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: editMaf.pl -i <file> -o <file> [OPTIONS]\n";
	print STDERR " -i <file>        [maf file | STDIN]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -r <string>      [reference genome (default: hg19)]\n";
	print STDERR " -l <int>         [compute properties of MAFs > length (default: largest MAF block)]\n";
	print STDERR " -a               [compute average property of all MAFs]\n";
	print STDERR " -s <int>         [compute secondary structure (0: RNAz, 1: CMfinder, 2: both)]\n";
	print STDERR " -h <help>\n";
	exit(-1);
}

##################################################################################################

## open coordinate file in BED format
my @data=();
if(defined($mafFile)) {
	my $size= -s "$mafFile";
	$size=$size/(1024*1024*1024);
	if($size>=1) {
		tie @data, 'Tie::File', $mafFile or die $!;
	}
	else {
		@data=openFile($mafFile);
	}
}
else {
	@data=openFile();
}

## determine maximum MAF block length in reference
## genome ($refGenome)
my $MAX_LENGTH=0;
foreach my $l(grep(/^s $refGenome\./, @data)) {
	my @t=split(/\s+/,$l);
	if($t[3]>$MAX_LENGTH) {
		$MAX_LENGTH=$t[3];
	}
}
## set length threshold to maximum length, if $threshold{'len'} is not defined
$threshold{'len'}=$MAX_LENGTH if(!defined($threshold{'len'}));
#print "$threshold{'len'}\n"; exit;

## read MAF block and compute properties
my %genome=(); my %mafProp=(); my $MAF_COUNTER=0;
my %mafInfo=(); my @comment=(); my $result=();
foreach my $l(@data) {
	chomp($l);
	if($l=~/^\#/) { push(@comment, $l); }
	elsif($l=~/^a score/) {
		if(defined($mafInfo{'aln'}) && $mafInfo{'len'}{$refGenome}>=$threshold{'len'}) {
			$MAF_COUNTER++;
			%{$mafProp{$MAF_COUNTER}}=compMafProp(%mafInfo);
			
			## print results, if $computeAvg==0
			if(!defined($computeAvg)) {
				if(defined($mafFile)) {
					printf("%s\t%0.2f\t%0.2f\t%0.2f\t", $mafFile, $mafProp{$MAF_COUNTER}{'avgPercentIdentity'}, $mafProp{$MAF_COUNTER}{'avgMAFLength'}, $mafProp{$MAF_COUNTER}{'avgGenomeCount'});
				}
				else {
					printf("STDIN\t%0.2f\t%0.2f\t%0.2f\t", $mafProp{$MAF_COUNTER}{'avgPercentIdentity'}, $mafProp{$MAF_COUNTER}{'avgMAFLength'}, $mafProp{$MAF_COUNTER}{'avgGenomeCount'});
				}

				## print secondary structure information
				if(defined($computeSS)) {
					if($computeSS==0 || $computeSS==2) {
						print "$mafProp{$MAF_COUNTER}{'sigValue'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'probability'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'consSeq'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'ss'}{'rnaz'}\t"
					}
					if($computeSS==1 || $computeSS==2) {
						print "$mafProp{$MAF_COUNTER}{'sigValue'}{'cmfinder'}\t$mafProp{$MAF_COUNTER}{'consSeq'}{'cmfinder'}\t$mafProp{$MAF_COUNTER}{'ss'}{'cmfinder'}\t";
					}
				}
				else { print "NA\tNA\tNA\tNA\tNA\tNA\tNA"; }
				print "\n";
			}
		}
		%mafInfo=();
		@{$mafInfo{'comment'}}=@comment;
		@comment=();
		$mafInfo{'score'}=$l;
	}
	elsif($l=~/^s\s+/) {
		my @t=split(/\s+/, $l);
		my($genome, $chr) = split(/\./, $t[1]);
		$mafInfo{'chr'}{$genome}=$chr;
		$mafInfo{'start'}{$genome}=$t[2];
		$mafInfo{'end'}{$genome}=$t[2]+($t[3]-1);
		$mafInfo{'len'}{$genome}=$t[3];
		$mafInfo{'strand'}{$genome}=$t[4];
		$mafInfo{'srcSize'}{$genome}=$t[5];
		$mafInfo{'aln'}{$genome}="$t[6]";
		push(@{$genome{'all'}}, $genome);
	}
}

## read last MAF block and compute properties
if(defined($mafInfo{'aln'}) && $mafInfo{'len'}{$refGenome}>=$threshold{'len'}) {
	$MAF_COUNTER++;
	%{$mafProp{$MAF_COUNTER}}=compMafProp(%mafInfo);
	%mafInfo=();

	## print results, if $computeAvg==0
	if(!defined($computeAvg)) {
		#printf("%s\t%0.2f\t%0.2f\t%0.2f\n", $mafProp{$MAF_COUNTER}{'mafId'}, $mafProp{$MAF_COUNTER}{'avgPercentIdentity'}, $mafProp{$MAF_COUNTER}{'avgMAFLength'}, $mafProp{$MAF_COUNTER}{'avgGenomeCount'});
		if(defined($mafFile)) {
			printf("%s\t%0.2f\t%0.2f\t%0.2f\t", $mafFile, $mafProp{$MAF_COUNTER}{'avgPercentIdentity'}, $mafProp{$MAF_COUNTER}{'avgMAFLength'}, $mafProp{$MAF_COUNTER}{'avgGenomeCount'});
		}
		else {
			printf("STDIN\t%0.2f\t%0.2f\t%0.2f\t", $mafProp{$MAF_COUNTER}{'avgPercentIdentity'}, $mafProp{$MAF_COUNTER}{'avgMAFLength'}, $mafProp{$MAF_COUNTER}{'avgGenomeCount'});
		}

		## print secondary structure information
		if(defined($computeSS)) {
			if($computeSS==0 || $computeSS==2) {
				print "$mafProp{$MAF_COUNTER}{'sigValue'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'probability'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'consSeq'}{'rnaz'}\t$mafProp{$MAF_COUNTER}{'ss'}{'rnaz'}\t"
			}
			if($computeSS==1 || $computeSS==2) {
				print "$mafProp{$MAF_COUNTER}{'sigValue'}{'cmfinder'}\t$mafProp{$MAF_COUNTER}{'consSeq'}{'cmfinder'}\t$mafProp{$MAF_COUNTER}{'ss'}{'cmfinder'}\t";
			}
		}
		else { print "NA\tNA\tNA\tNA\tNA\tNA\tNA"; }
		print "\n";
	}
}

## print results, if compute average of MAF properties
if($computeAvg && $MAF_COUNTER>0) {
	if(defined($mafFile)) {	print "$mafFile\t"; }
	else { print "STDIN\t"; }

	my $percentIdentity=0; my $mafLength=0; my $genomeCount=0;
	foreach my $maf(keys(%mafProp)) {
		## compute average percent identity
		$percentIdentity+=$mafProp{$maf}{'avgPercentIdentity'};

		## compute average maf length
		$mafLength+=$mafProp{$maf}{'avgMAFLength'};

		## compute average genome count
		$genomeCount+=$mafProp{$maf}{'avgGenomeCount'};
	}

	printf("%0.2f\t%0.2f\t%0.2f\n", $percentIdentity/$MAF_COUNTER, $mafLength/$MAF_COUNTER, $genomeCount/$MAF_COUNTER);
	print $result if(defined($computeSS));
}
elsif($MAF_COUNTER==0) {
	if(defined($mafFile)) {	print "$mafFile\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"; }
	else { print "STDIN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"; }
}

exit(0);

#######################################################################
sub compMafProp {
	my($mafInfo) = @_;

	my %mafProp=();
	## open temporary file to keep fasta sequence per MAF block
	my $tmpFileName=$refGenome."_".$mafInfo{'chr'}->{$refGenome}."_".$mafInfo{'start'}->{$refGenome}."_".$mafInfo{'len'}->{$refGenome};
	open(TMPFILE, ">$tmpFileName.fasta") || die $!;
	#print "$tmpFileName\n";

	## MAF identifier
	$mafProp{'mafId'}=$tmpFileName;

	## MAF block length
	$mafProp{'avgMAFLength'}=$mafInfo{'len'}{$refGenome};

	foreach my $gen(keys(%{$mafInfo{'aln'}})) {
		if(defined($computeSS) && $computeSS==1) {
			$mafInfo{'aln'}->{$gen}=~s/\-//g;
			my @t=split(//,$mafInfo{'aln'}->{$gen});
			if(scalar(@t)>40) {
				print TMPFILE ">$gen"."_".$mafInfo{'chr'}->{$gen}."_".$mafInfo{'start'}->{$gen}."_".$mafInfo{'len'}->{$gen}."\n";
				print TMPFILE "$mafInfo{'aln'}->{$gen}\n";
			}
		}
		else {
			print TMPFILE ">$gen"."_".$mafInfo{'chr'}->{$gen}."_".$mafInfo{'start'}->{$gen}."_".$mafInfo{'len'}->{$gen}."\n";
			print TMPFILE "$mafInfo{'aln'}->{$gen}\n";
		}
	}
	close(TMPFILE);

	## compute average percent sequence identity in MAF block
	$mafProp{'avgPercentIdentity'}=`seqIdentity.pl -f $tmpFileName.fasta | tail -n 1`;
	chomp($mafProp{'avgPercentIdentity'});
	#print "$mafProp{'avgPercentIdentity'}\n";

	## number of genomes included in the MAF block
	$mafProp{'avgGenomeCount'}=scalar(keys(%{$mafInfo{'aln'}}));

	## Compute secondary structure prediction using
	## RNAz
	if(defined($computeSS) && ($computeSS==0 || $computeSS==2)) {
		## PETCoFold
		#$ENV{PETFOLDBIN}='/home/users/seemann/projects/PETfold/code/c/bin/';
		#push(@SS, `PETfold -f $tmpFileName.fasta -g 0.50`);

		system("fasta2clustal.pl -i $tmpFileName.fasta -o $tmpFileName.clustalw -f clustalw >/dev/null");

		## RNAAliFold
		#$mafProp{'SS'}=`RNAalifold -p -r -d2 -noLP -color -aln < $tmpFileName.clustalw`;

		my @data=`RNAz $tmpFileName.clustalw 2>/dev/null`;
		my $start=0;
		foreach(@data) {
			chomp($_);
			if($_=~/Mean z\-score\:/) {
				$_=~s/^.+score\:\s+//g;
				$mafProp{'sigValue'}{'rnaz'}=$_;
			}
			elsif($_=~/SVM RNA\-class probability/) {
				$_=~s/^.+probability\:\s+//g;
				$mafProp{'probability'}{'rnaz'}=$_;
			}
			elsif($start==1 && $_=~/^$/) {
				$start=0;
			}
			elsif($_=~/\>consensus/) {
				$start=1;
			}
			elsif($start==1 && $_=~/^[A-Z\-\_]+/i) {
				$mafProp{'consSeq'}{'rnaz'}.=$_;
			}
			elsif($start==1 && $_=~/\(/) {
				$_=~s/\s+.+//g;
				$mafProp{'ss'}{'rnaz'}.=$_;
			}
		}

		$mafProp{'sigValue'}{'rnaz'}="NA" if(!defined($mafProp{'sigValue'}{'rnaz'}));
		$mafProp{'consSeq'}{'rnaz'}="NA" if(!defined($mafProp{'consSeq'}{'rnaz'}));
		$mafProp{'ss'}{'rnaz'}="NA" if(!defined($mafProp{'ss'}{'rnaz'}));

		## remove temporary file
		system("rm $tmpFileName.clustalw");
		#system("rm alifold.out alidot.ps aln.ps");
		#system("mv alirna.ps $tmpFileName.ps");

		#print "$mafProp{'sigValue'}{'rnaz'}\t$mafProp{'consSeq'}{'rnaz'}\t$mafProp{'ss'}{'rnaz'}\n";
	}
	else {
		$mafProp{'sigValue'}{'rnaz'}="NA";
		$mafProp{'consSeq'}{'rnaz'}="NA";
		$mafProp{'ss'}{'rnaz'}="NA";
	}

	## CMfinder
	if(defined($computeSS) && ($computeSS==1 || $computeSS==2)) {
		if(! -e $tmpFileName) {
			$ENV{CMfinder}='/usr/local/opt/cmfinder/CMfinder_0.2.2/bin';
			$ENV{BLAST}='/usr/local/bin';
			system("\$CMfinder/cmfinder.pl $tmpFileName.fasta >/dev/null");
			system("mkdir $tmpFileName");
			system("mv $tmpFileName.fasta* $tmpFileName");
			system("mv latest.cm $tmpFileName");
			$ENV{LD_LIBRARY_PATH}='/chromo/4/projects/users/sdu/elfar/software/CMfinder_03/lib';
			$ENV{CMfinder}='/chromo/4/projects/users/sdu/elfar/software/CMfinder_03/';
			system("for i in $tmpFileName/*.motif.*; do \$CMfinder/bin/posterior --partition \$i > \$i.pscore; done");
		}

		## parse output
		if(glob("$tmpFileName/*pscore")) {
			my @data=`cat $tmpFileName/*pscore`;
			my %result=(); my $i=0;
			foreach(@data) {
				chomp($_);
				if($_=~/SS\_cons/) {
					$_=~s/^.+SS\_cons\s+//g;
					$result{$i}{'ss'}.=$_;
				}
				elsif($_=~/^$refGenome\_/) {
					$_=~s/^$refGenome.+\s+//g;
					$result{$i}{'consSeq'}.=$_;
				}
				elsif($_=~/Total pair posterior/) {
					$_=~s/Total pair posterior\s+//g;
					$result{$i}{'pscore'}=$_;
					$i++;
				}
			}

			foreach(reverse sort { $result{$a}{'pscore'} <=> $result{$b}{'pscore'} } keys(%result)) {
				$mafProp{'ss'}{'cmfinder'}=$result{$_}{'ss'};
				$mafProp{'consSeq'}{'cmfinder'}=$result{$_}{'consSeq'};
				$mafProp{'sigValue'}{'cmfinder'}=$result{$_}{'pscore'};
				last;
			}
		}
		
		$mafProp{'ss'}{'cmfinder'}="NA"	if(!defined($mafProp{'ss'}{'cmfinder'}));
		$mafProp{'consSeq'}{'cmfinder'}="NA" if(!defined($mafProp{'consSeq'}{'cmfinder'}));
		$mafProp{'sigValue'}{'cmfinder'}="NA" if(!defined($mafProp{'sigValue'}{'cmfinder'}));

		#print "$mafProp{'sigValue'}{'cmfinder'}\t$mafProp{'consSeq'}{'cmfinder'}\t$mafProp{'ss'}{'cmfinder'}\n";
	}
	else {
		$mafProp{'ss'}{'cmfinder'}="NA";
		$mafProp{'consSeq'}{'cmfinder'}="NA";
		$mafProp{'sigValue'}{'cmfinder'}="NA";
	}

	system("rm $tmpFileName.fasta");

	#print "$mafProp{'avgPercentIdentity'}\t$mafProp{'avgMAFLength'}\t$mafProp{'avgGenomeCount'}\n";
	return %mafProp;
}

## convert maf to fasta
sub maf2fasta {
	my $tmpFileName=$refGenome."_".$mafInfo{'chr'}->{$refGenome}."_".$mafInfo{'start'}->{$refGenome}."_".$mafInfo{'len'}->{$refGenome};
	open(TMPFILE, ">$tmpFileName.fasta") || die $!;
	#print "$tmpFileName\n";

	foreach my $gen(keys(%{$mafInfo{'aln'}})) {
		if(defined($computeSS) && $computeSS==1) {
			$mafInfo{'aln'}->{$gen}=~s/\-//g;
			my @t=split(//,$mafInfo{'aln'}->{$gen});
			if(scalar(@t)>40) {
				print TMPFILE ">$gen"."_".$mafInfo{'chr'}->{$gen}."_".$mafInfo{'start'}->{$gen}."_".$mafInfo{'len'}->{$gen}."\n";
				print TMPFILE "$mafInfo{'aln'}->{$gen}\n";
			}
		}
		else {
			print TMPFILE ">$gen"."_".$mafInfo{'chr'}->{$gen}."_".$mafInfo{'start'}->{$gen}."_".$mafInfo{'len'}->{$gen}."\n";
			print TMPFILE "$mafInfo{'aln'}->{$gen}\n";
		}
	}
	close(TMPFILE);
}

## stitch MAF blocks
sub stitchMaf {
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
