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

use vars qw ($help $bgFile $annoFile $overlap $sizeCutoff $printOpt $progDir $flag);

##########################################################################################################
## parse input options
$overlap = 0.000000000000001;
$sizeCutoff = -1;
$printOpt = 0;
$flag=0;

GetOptions ("c=s"  => \$bgFile,
            "a=s"  => \$annoFile,
            "o=s"  => \$overlap,
            "s=s"  => \$sizeCutoff,
            "p=s"  => \$printOpt,
            "r=s"  => \$progDir,
            "f=s"  => \$flag,
            "help" => \$help,
            "h"    => \$help);

usage() if ($help || !$bgFile || !$annoFile);

##########################################################################################################

sub usage {
	print STDERR "\nusage: flagCluster.pl -c <file> -a <file> [OPTIONS]\n";
	print STDERR "flag block groups overlapping with annotated ncRNAs or loci\n";
	print STDERR "\n";
	print STDERR "[INPUT]\n";
	print STDERR " -c <file>    [block group file]\n";
	print STDERR " -a <file>    [ncRNAs or loci annotation file(s) seperated by comma or directory]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <file>    [min overlap in % (default = 0)]\n";
	print STDERR " -s <file>    [size cutoff, i.e. 1.5 times of the length (default = off)]\n";
	print STDERR " -p <int>     [print option (0 = all [default]; 1 = only block groups with annotation)]\n";
	print STDERR " -r <dir>     [directory with programs (coor2annotation.pl)]\n";
	print STDERR " -f <int>     [flag (0 = annotation [default]; 1 = loci)]\n";
	print STDERR "[VERSION]\n";
	print STDERR " 15-01-2014\n";
	print STDERR "[BUGS]\n";
	print STDERR " Please report bugs to sachin\@rth.dk\n";
	print STDERR "\n";
	exit(-1);
}

##########################################################################################################

## retrieve annotation file(s) and path
my %file=(); my $inputIsDir=();

if(-d $annoFile) {
	$inputIsDir=1;
	if($flag==0) {
		@{$file{'name'}}=`ls $annoFile | perl -ane 'chomp(\$_); \$_=~s/format\..+\$/format/g; if(\$_=~/ncRNAs/ || \$_=~/rfam/) { \$seen{\$_}=1; } END { foreach(keys(%seen)) { print "\$_\n"; } }'`;
	}
	elsif($flag==1) {
		@{$file{'name'}}=`ls $annoFile | perl -ane 'chomp(\$_); \$_=~s/format\..+\$/format/g; if(\$_=~/loci/) { \$seen{\$_}=1; } END { foreach(keys(%seen)) { print "\$_\n"; } }'`;
	}
}
else {
	$annoFile=~m/^.+\//g;
	$file{'path'}=$&;
	if(!defined($file{'path'})) { $file{'path'}="."; }
	foreach(split(/\,/, $annoFile)) {
		$_=~s/^.+\///g;
		if(-e "$file{'path'}/$_") {
			push(@{$file{'name'}}, $_);
		}
	}
}

## open block group file
if($bgFile=~/\.gz$/) { open(FILE, "gunzip -c $bgFile |") || die $!; }
else { open(FILE, "<$bgFile") || die "cannot open $bgFile\n"; }

## print block group 
if(!defined($file{'name'}) && ! -d $annoFile) {
	while(<FILE>) {
		chomp;
		my @t=split(/\s+/, $_);
		if($_=~/^\>/ && scalar(@t)==8) {
			print "$_\tn/a\tn/a\tn/a\n";
		}
		else { print "$_\n"; }
	}
	exit(0);
}

if($flag==0) {
	## open annotation file
	my %data=();
	if(defined($inputIsDir)) {
		my @temp=`zless $annoFile/* | sort -k 2n,2 -k 3n,3`;
		foreach my $l(@temp) {
			my @t=split(/\s+/, $l);
			push(@{$data{$t[0]}}, $l);
		}
	}
	else {
		my $file=();
		foreach(@{$file{'name'}}) {
			chomp($_);
			$file.="$file{'path'}/$_ ";
		}
		my @temp=`zless $file | sort -k 2n,2 -k 3n,3`;
		foreach my $l(@temp) {
			my @t=split(/\s+/, $l);
			push(@{$data{$t[0]}}, $l);
		}
	}

	## exit, if annotation file is empty
	if(keys(%data)==0) {
		print STDERR "Input annotation file is empty.\n";
		exit(-1);
	}

	## retrieve annotations overlapping to the coordinates of the block groups
	my $s = 0;
	while(<FILE>){
		chomp;
		next if(/^\#/);
		my @line = split(/\s+/, $_); my $doFlag=0;
		if(/^\>/ && scalar(@line)==8) { $doFlag=1; }
		elsif(/^\>/ && scalar(@line)==11 && $line[9]=~/n\/a/) { $doFlag=1; }
		elsif(/^\>/ && scalar(@line)!=11 && scalar(@line)!=8) { print STDERR "Incorrect input format of $bgFile"; exit; }

		if(/^\>/ && $doFlag){
			my $flag = "a"; $s = 0;
			my $perOverlap=0;
			my $anno=();
			if(defined($progDir)) {
				if(defined($inputIsDir)) {
					foreach(@{$file{'name'}}) {
						chomp($_);
	 					if($line[4]=~/^\+$/ && -e "$annoFile/$_.$line[1].P.gz") {
							$anno=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].P.gz -s $line[4] -o $overlap`;
						}
						elsif($line[4]=~/^\-$/ && -e "$annoFile/$_.$line[1].N.gz") {
							$anno=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].N.gz -s $line[4] -o $overlap`
						}
						else { $anno="No"; }
						last if($anno!~/No/);
					}
				}
				else {
					foreach(@{$file{'name'}}) {
						chomp($_);
						$anno=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $file{'path'}/$_ -s $line[4] -o $overlap`;
						last if($anno!~/No/);
					}
				}
			}
			else {
				$anno=coor2annotation("$line[1]:$line[2]-$line[3]", $line[4], \@{$data{$line[1]}}, $overlap);
=cu
				if(defined($inputIsDir)) {
					foreach(@{$file{'name'}}) {
						chomp($_);
						if($line[4]=~/^\+$/ && -e "$annoFile/$_.$line[1].P.gz") {
							$anno=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].P.gz -s $line[4] -o $overlap`;
						}
						elsif($line[4]=~/^\-$/ && -e "$annoFile/$_.$line[1].N.gz") {
							$anno=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].N.gz -s $line[4] -o $overlap`;
						}
						else { $anno="No"; }
						last if($anno!~/No/);
					}
				}
				else {
					foreach(@{$file{'name'}}) {
						chomp($_);
						$anno=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $file{'path'}/$_ -s $line[4] -o $overlap`;
						last if($anno!~/No/);
					}
				}
=cut
			}

			if(defined($anno) && $anno!~/^No/) {
				my @t=split(/\s+/, $anno);
				$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$t[3]\t$t[7]\t$t[8]";
				$s = 1;
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
	## open loci file
	my %data=();
	if(defined($inputIsDir)) {
		my @temp=`zless $annoFile/* | sort -k 2n,2 -k 3n,3`;
		foreach my $l(@temp) {
			my @t=split(/\s+/, $l);
			push(@{$data{$t[0]}}, $l);
		}
	}
	else {
		my $file=();
		foreach(@{$file{'name'}}) {
			chomp($_);
			$file.="$file{'path'}/$_ ";
		}
		my @temp=`zless $file | sort -k 2n,2 -k 3n,3`;
		foreach my $l(@temp) {
			my @t=split(/\s+/, $l);
			push(@{$data{$t[0]}}, $l);
		}
	}

	## exit, if annotation file is empty
	if(keys(%data)==0) {
		print STDERR "Input annotation file is empty.\n";
		exit(-1);
	}

	my $s = 0;
	while(<FILE>){
		chomp;
		next if(/^\#/);
		if(/^\>/){
			my @line = split(/\s+/, $_);
			if(scalar(@line)!=11) { print STDERR "Input file is not flagged with ncRNA annotation\nLine: $_\n"; exit; }
			my $flag = "a"; $s = 0;
			my $perOverlap=0;
			my $loci=();
			if(defined($progDir)) {
				if(defined($inputIsDir)) {
					foreach(@{$file{'name'}}) {
						chomp($_);
 						if($line[4]=~/^\+$/ && -e "$annoFile/$_.$line[1].P.gz") {
							$loci=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].P.gz -s $line[4] -o $overlap`;
						}
						elsif($line[4]=~/^\-$/ && -e "$annoFile/$_.$line[1].N.gz") {
							$loci=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].N.gz -s $line[4] -o $overlap`
						}
						else { $loci="No"; }
						last if($loci!~/^No/);
					}
				}
				else {
					foreach(@{$file{'name'}}) {
						chomp($_);
						$loci=`$progDir/coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $file{'path'}/$_ -s $line[4] -o $overlap`;
						last if($loci!~/^No/);
					}
				}
			}
			else {
				$loci=coor2annotation("$line[1]:$line[2]-$line[3]", $line[4], \@{$data{$line[1]}}, $overlap);
=cu
				if(defined($inputIsDir)) {
					foreach(@{$file{'name'}}) {
						chomp($_);
 						if($line[4]=~/^\+$/ && -e "$annoFile/$_.$line[1].P.gz") {
							$loci=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].P.gz -s $line[4] -o $overlap`;
						}
						elsif($line[4]=~/^\-$/ && -e "$annoFile/$_.$line[1].N.gz") {
							$loci=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $annoFile/$_.$line[1].N.gz -s $line[4] -o $overlap`
						}
						else { $loci="No"; }
						last if($loci!~/^No/);
					}
				}
				else {
					foreach(@{$file{'name'}}) {
						chomp($_);
						$loci=`coor2annotation.pl -c \'$line[1]:$line[2]-$line[3]\' -f $file{'path'}/$_ -s $line[4] -o $overlap`;
						last if($loci!~/^No/);
					}
				}
=cut
			}
			if(defined($loci) && $loci!~/^No/ && ($line[10]=~/\./ || $line[10]=~/n\/a/)) {
				my @t=split(/\s+/, $loci);
				$flag = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$t[7]";
				$s = 1;
			}
			elsif(defined($loci) && $loci!~/^No/) {
				my @t=split(/\s+/, $loci);
				$line[10] = $t[7]."_".$line[10];
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

## retrieve non-coding RNA annotation overlapping to the input coordinate (sorted by start coordinate in BED format)
sub coor2annotation {
	my %coor=();
	my $annotation=();
	my $minPerOverlap=();
	($coor{'coor'}, $coor{'strand'}, $annotation, $minPerOverlap)=@_;

	($coor{'chr'}, $coor{'start'}, $coor{'end'}) = split(/[\:\-]+/, $coor{'coor'});

	if(!defined($coor{'chr'}) || !defined($coor{'start'}) || !defined($coor{'end'})) { print STDERR "Incorrect coordinate\n"; return "None"; }

	if(scalar(@{$annotation})==0) {
		return "None";
	}

	my $posmin=0;
	#my $posmax=$#annotation;
	my $posmax=scalar(@{$annotation})-1;
	my $result=0;
	while($result==0) {
		my $mid = int (($posmin + $posmax) /2);
		#print "$annotation->[$posmin]\n$annotation->[$mid]\n$annotation->[$posmax]\n\n";
		$result = binarySearchCoor($coor{'start'}, $coor{'end'}, $mid, \$posmin, \$posmax, "bed", $annotation->[$mid]);
		$result=1 if($posmin==$posmax);
	}

	$coor{'annoOverlap'}=$minPerOverlap;
	$coor{'annoLength'}=10000000000000000000000000000000000000000000000;

	#foreach my $l(@{$annotation[$posmin..$posmax]}) {
	for(my $i=$posmin; $i<=$posmax; $i++) {
		my @anno=split(/\s+/, $annotation->[$i]);
		my $overlap=(); my $length=();

		if(defined($coor{'strand'})) {
			$overlap=checkOverlap($anno[1], $anno[2], $coor{'start'}, $coor{'end'}, 1, $anno[0], $coor{'chr'}, $anno[5], $coor{'strand'});
			#print "$anno[1]\t$anno[2]\t$coor{'start'}\t$coor{'end'}\t$anno[0]\t$coor{'chr'}\t$anno[5]\t$coor{'strand'}\t$anno[3]\t$overlap\t";
		}
		else {
			$overlap=checkOverlap($anno[1], $anno[2], $coor{'start'}, $coor{'end'}, 1, $anno[0], $coor{'chr'});
		}
		$length=($anno[2]-$anno[1])+1;
		#print "$length\n";

		if($overlap >= $coor{'annoOverlap'} && $length < $coor{'annoLength'}) {
			$coor{'anno'}=$annotation->[$i];
			$coor{'annoOverlap'}=$overlap;
			## determine the length of the overlapping annotation
			$coor{'annoLength'}=($anno[2]-$anno[1])+1;
			#print "$coor{'anno'}\n";
		}
	}

	## print the annotation overlapping to the input coordinate
	if(defined($coor{'anno'})) {
		return "$coor{'anno'}";
	}
	else {
		return "None";
	}
}

## check overlap between two coordinates
sub checkOverlap {
	my ($tStart, $tEnd, $bStart, $bEnd, $percentage, $tChr, $bChr, $tStrand, $bStrand) = @_;
	return 0 if(!defined($tStart) || !defined($tEnd) || !defined($bStart) || !defined($bEnd));
	if(defined($tChr) && defined($bChr)) {
		return 0 if($tChr!~/^$bChr$/i);
	}
	if(defined($tStrand) && defined($bStrand)) {
		return 0 if($tStrand!~/\Q$bStrand\E/);
	}

	my $overlap=0;
	for(my $k=$tStart; $k<=$tEnd; $k++) {
		if($k>=$bStart && $k<=$bEnd) {
			$overlap++;
		}
	}

	## take first coordinate as reference
	if(defined($percentage) && $percentage==2) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	## take second coordinate as reference
	elsif(defined($percentage) && $percentage==3) {
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}
	## take shortest coordinate as reference
	elsif(defined($percentage) && $percentage==1 && $tEnd-$tStart <= $bEnd-$bStart) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	elsif(defined($percentage) && $percentage==1) {
		if((($bEnd-$bStart)+1)==0) { print STDERR "Error, start=$bStart and end=$bEnd\n"; exit; }
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}
	## return 1 if overlap, return 0 otherwise
	elsif(defined($percentage) && $percentage==4 && $overlap>0) {
		$overlap=1;
	}
	#print "$tStart\t$tEnd\t$bStart\t$bEnd\t$overlap\n";
	return $overlap;
}

## determine overlapping coordinates by binary sort (coordinate file should be sorted)
sub binarySearchCoor {
	my($start, $end, $mid, $posmin, $posmax, $mapFormat, $line)=@_;
	my $result=();
	my @t=split(/\s+/, $line);
	if($mapFormat=~/segemehl/) {
		if($t[10]>$end) { ${$posmax}=$mid; $result=0; }
		elsif($t[11]<$start) { ${$posmin}=$mid; $result=0; }
		else { $result=1; }
	}
	elsif($mapFormat=~/bed/) {
		if($t[1]>$end) { ${$posmax}=$mid; $result=0; }
		elsif($t[2]<$start) { ${$posmin}=$mid; $result=0; }
		else { $result=1; }
	}

	if((${$posmax}-${$posmin})==1) { $result=1; }
	return $result;
}

