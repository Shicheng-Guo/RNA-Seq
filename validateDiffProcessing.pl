#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;
use List::Compare;
use Data::Dumper;
use Tie::IxHash;

#######################################################################################################
use vars qw($compute $tissues $dpLociFile $exprLociFile $overlapFile $bedFile $win $steps $deseqMatrixFile $bgFile $bins $mapDir $treeDir $sizeFactorFile $ucscDir $webPath $replicates $rep1UCSCDir $rep2UCSCDir $overlapFileRep1Rep2 $overlapFileRep1All $dbaScoreDiffTissue $bgExprLociFile $bgDir $removeDummy $refGenome $outDir $id $help);
## GLOBAL VARIABLES
my %threshold=();
$threshold{'minUniqReads'}=0;
$threshold{'minClusterScore'}=0.15;
$threshold{'minFScore'}=1;
$threshold{'minTissuesInCluster'}=2;
$threshold{'consistentClusters'}=6;
my $minBlockOverlap=90;
my $mapFormat="segemehl";
$webPath="http://rth.dk/~sachin/ucsc_tracks/";
$refGenome="/home/projects/deepalign/data/genomes/hg19/hg19.fa";
$outDir=".";

GetOptions ("c=s"  => \$compute,
            "z=i"  => \$tissues,
            "f=s"  => \$dpLociFile,
            "s=s"  => \$exprLociFile,
            "r=i"  => \$replicates,
            "o=s"  => \$overlapFile,
            "e=s"  => \$bedFile,
            "w=i"  => \$win,
            "p=i"  => \$steps,
            "d=s"  => \$deseqMatrixFile,
            "k=s"  => \$bgFile,
            "b=i"  => \$bins,
            "m=s"  => \$mapDir,
            "t=s"  => \$treeDir,
            "v=s"  => \$sizeFactorFile,
            "u=s"  => \$ucscDir,
            "a=s"  => \$webPath,
            "i=s"  => \$threshold{'minClusterScore'},
            "j=i"  => \$threshold{'minFScore'},
            "n=i"  => \$threshold{'minTissuesInCluster'},
            "q=i"  => \$threshold{'minUniqReads'},
            "o1=s" => \$overlapFileRep1Rep2,
            "o2=s" => \$overlapFileRep1All,
            "g=i"  => \$threshold{'consistentClusters'},
            "l"    => \$dbaScoreDiffTissue,
            "x=s"  => \$bgExprLociFile,
            "y=s"  => \$bgDir,
            "y1"   => \$removeDummy,
            "g1"   => \$refGenome,
            "o3=s" => \$outDir,
						"o4=s" => \$id,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$compute || !$tissues);

#######################################################################################################
sub usage {
	print STDERR "\nProgram: validateDiffProcessing.pl (validate differential processing)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: validateDiffProcessing.pl -c <computation option> [OPTIONS]\n";
	print STDERR "          -c <computation option>\n";
	print STDERR "          -z <number of tissues analyzed>\n";
	print STDERR "          -h <help>\n";
	print STDERR "Option a: convert to BED format and/or determine overlap between replicates\n";
	print STDERR "          -s <file>         [file with a list of consitutively expressed loci]\n";
	print STDERR "          -q <int>          [print loci with reads > minimum unique reads (default: 0)]\n";
	print STDERR "          -r <int>          [perform analysis for replicates (1: rep1 with 2; 2: rep1 or 2 with all]\n";
	print STDERR "          -o <file>         [output file with overlapping coordinates and cluster information]\n";
	print STDERR "Option b: compute sequence bias (mono, di and tri) in nucleotide sequence(s)\n";
	print STDERR "          -e <file>         [BED file | STDIN]\n";
	print STDERR "          -w <int>          [window]\n";
	print STDERR "          -p <int>          [steps]\n";
	print STDERR "Option c: check for differential expression in consitutively expressed loci\n";
	print STDERR "          -f <file>         [file with a list of differentially expressed loci]\n";
	print STDERR "          -d <file>         [DESeq matrix file]\n";
	print STDERR "Option d: read position bias in block group(s)\n";
	print STDERR "          -k <file>         [block group file in blockbuster format | STDIN]\n";
	print STDERR "Option e: position bias in block group(s) relative to genes\n";
	print STDERR "          <STDIN>           [intersection of block group and gene coordinates in BED format]\n";
	print STDERR "                            [FORMAT: chr start end name score strand chr start end name score strand]\n";
	print STDERR "          -b <int>          [bins]\n";
	print STDERR "Option f: create and visualize UCSC tracks for differentially processed loci\n";
	print STDERR "          -o <file>         [input file with overlapping coordinate in BED format]\n";
	print STDERR "          -t <dir>          [directory with block group, alignment and tree files for all loci]\n";
	print STDERR "          -v <file>         [input file having size factors]\n";
	print STDERR "          -y1               [do not plot dummy blocks]\n";
	print STDERR "          -o4 <string>      [run for specific id]\n";
	print STDERR "Option g: identify a list of validated differentially processed loci accross replicates (if present)\n";
	print STDERR "          -f <file>         [file with a list of differentially processed loci]\n";
	print STDERR "          -o <file>         [file with overlapping coordinates (-r), or;]\n";
	print STDERR "                            [constitutively expressed loci file in BED format]\n";
	print STDERR "          -i <int>          [print loci with >= minimum cluster score (default: 0.15)]\n";
	print STDERR "          -j <int>          [print loci with >= minimum fScore (default: 1)]\n";
	print STDERR "          -n <int>          [print loci with >= minimum tissues in a cluster (default: 2)]\n";
	print STDERR "                            [recommended: 4 (both replicates clustered together) or 2 (otherwise)]\n";
	print STDERR "          -r <int>          [perform analysis for replicates]\n";
	print STDERR "                            [1: rep1 with 2]\n";
	print STDERR "                            [2: rep1 or 2 with all]\n";
	print STDERR "                            [3: only all]\n";
	print STDERR "                            [4: only rep1 or rep2]\n";
	print STDERR "Option h: correlate DBA score and consistent clustering between replicates\n";
	print STDERR "          -o1 <file|STDIN>  [overlap file between replicate 1 and 2 (eg: rep1_rep2.overlap.valid)\n";
	print STDERR "          -o2 <file>        [overlap file between replicate 1 and all (default: repl1_all.overlap)\n";
	print STDERR "          -g <int>          [cut-off for consistent clustering (default: 6)]\n";
	print STDERR "          -t <dir>          [directory with block group, alignment and tree files for all loci]\n";
	print STDERR "          -l                [print DBA score between profiles from different tissues]\n";
	print STDERR "Option i: convert bgExpCons to BED format\n";
	print STDERR "          -x                [file with block group information at each loci]\n";
	print STDERR "          -y                [block group directory]\n";
	print STDERR "Option j: check for read position bias using seqbias package (R)\n";
	print STDERR "          -f <file>         [file with a list of differentially expressed loci]\n";
	print STDERR "          -t <dir>          [directory with block group, alignment and tree files for all loci]\n";
	print STDERR "          -m <dir>          [input directory having read mapping files and yaml (seqbias model) files]\n";
	print STDERR "          -o3 <dir>         [output directory having normalized read counts (default: current)]\n";
	print STDERR "          -g1 <dir>         [reference genome sequence (default: /home/projects/deepalign/data/genomes/hg19/hg19.fa]\n\n";
	exit(-1);
}

#######################################################################################################

if($compute=~/[aA]+/) {
	usage() if(!$exprLociFile);

	## convert consitutively expressed loci file to BED, determine overlap between replicates
	if(defined($exprLociFile) && defined($replicates) && $replicates==1) {
		usage() if(!$overlapFile);

		## determine files corresponding to both replicate 1 and 2
		my %file=(); my %name=();
		if($exprLociFile=~/Rep1/) {
			$file{'rep1'}="$exprLociFile";
			$file{'rep2'}="$exprLociFile";
			$file{'rep2'}=~s/Rep1/Rep2/;
			if(! -e $file{'rep1'} || ! -e $file{'rep2'}) {
				print STDERR "Either of $file{'rep1'} or $file{'rep2'} do not exist\n";
				exit(-1);
			}
			$name{'rep1'}=$file{'rep1'}; $name{'rep1'}=~s/\.gz//g;
			$name{'rep2'}=$file{'rep2'}; $name{'rep2'}=~s/\.gz//g;
		}
		elsif($exprLociFile=~/Rep2/) {
			$file{'rep2'}="$exprLociFile";
			$file{'rep1'}="$exprLociFile";
			$file{'rep1'}=~s/Rep2/Rep1/;
			if(! -e $file{'rep1'} || ! -e $file{'rep2'}) {
				print STDERR "Either of $file{'rep1'} or $file{'rep2'} do not exist\n";
				exit(-1);
			}
			$name{'rep1'}=$file{'rep1'}; $name{'rep1'}=~s/\.gz//g;
			$name{'rep2'}=$file{'rep2'}; $name{'rep2'}=~s/\.gz//g;
		}
		else { print STDERR "Cannot locate both replicates for $exprLociFile\n"; exit(-1); }

		## check format of first replicate file
		my $found=`zless $file{'rep1'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=(14+$tissues)) { print STDERR "Incorrect format for $file{'rep1'}\n"; exit(-1); }

		## check format of second replicate file
		$found=`zless $file{'rep2'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=(14+$tissues)) { print STDERR "Incorrect format for $file{'rep2'}\n"; exit(-1); }

		## convert constitutively expressed loci file into bed format (replicate 1)
		system("echo \"#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\" > $name{'rep1'}.bed");
		system("zgrep -v \"\^\\#\" $file{'rep1'} | perl -an -F'/\\s+/' -e 'if(\$F[0]=~/bgP/) { \$strand=\"\+\"; } else { \$strand=\"\-\"; } \@t=split(/[\\:\\-]+/, \$F[1]); \$blocks=1; \$clusterInfo=(); \$uniqReads=100; foreach \$l(\@F[3..($tissues+2)]) { \@b=split(/\\|/,\$l); if(\$b[9]<\$uniqReads) { \$uniqReads=\$b[9]; } if(\$b[7]>\$blocks) { \$blocks=\$b[7]; } \$clusterInfo.=\"\$b[10],\" } \$clusterInfo=~s/\,\$//g; if(\$uniqReads>=$threshold{'minUniqReads'}) { print \"\$t[0]\\t\$t[1]\\t\$t[2]\\t\$F[0]\\t0\\t\$strand\\t\$F[$tissues+12]_\$F[$tissues+13]\\t\$F[$tissues+3]\\t\$F[$tissues+4]\\t\$F[$tissues+5]\\t\$F[$tissues+10]\\t\$F[2]\\t\$blocks\\t\$uniqReads\\t\$clusterInfo\\n\"; }' >> $name{'rep1'}.bed");

		## convert constitutively expressed loci file into bed format (replicate 2)
		system("echo \"#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\" > $name{'rep2'}.bed");
		system("zgrep -v \"\^\\#\" $file{'rep2'} | perl -an -F'/\\s+/' -e 'if(\$F[0]=~/bgP/) { \$strand=\"\+\"; } else { \$strand=\"\-\"; } \@t=split(/[\\:\\-]+/, \$F[1]); \$blocks=1; \$clusterInfo=(); \$uniqReads=100; foreach \$l(\@F[3..($tissues+2)]) { \@b=split(/\\|/,\$l); if(\$b[9]<\$uniqReads) { \$uniqReads=\$b[9]; } if(\$b[7]>\$blocks) { \$blocks=\$b[7]; } \$clusterInfo.=\"\$b[10],\" } \$clusterInfo=~s/\,\$//g; if(\$uniqReads>=$threshold{'minUniqReads'}) { print \"\$t[0]\\t\$t[1]\\t\$t[2]\\t\$F[0]\\t0\\t\$strand\\t\$F[$tissues+12]_\$F[$tissues+13]\\t\$F[$tissues+3]\\t\$F[$tissues+4]\\t\$F[$tissues+5]\\t\$F[$tissues+10]\\t\$F[2]\\t\$blocks\\t\$uniqReads\\t\$clusterInfo\\n\"; }' >> $name{'rep2'}.bed");

		## determine overlapping coordinates between the two replicates
		my @data=();
		if($exprLociFile=~/Rep1/) {
			#@data = `bed2intersect.pl -f $name{'rep1'}.bed -s $name{'rep2'}.bed | grep -v "^NA"`;
			@data = `intersectBed -a $name{'rep1'}.bed -b $name{'rep2'}.bed -wo -f 0.80 -r -s | cut -f 1-30`;
		}
		elsif($exprLociFile=~/Rep2/) {
			#@data = `bed2intersect.pl -f $name{'rep2'}.bed -s $name{'rep1'}.bed | grep -v "^NA"`;
			@data = `intersectBed -a $name{'rep2'}.bed -b $name{'rep1'}.bed -wo -f 0.80 -r -s | cut -f 1-30`;
		}

		## determine overlap for 'pvclust' clusters between the two replicates
		open(OUTFILE,">$overlapFile");
		print OUTFILE "#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\tchr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\tclusterOverlap\n";
		foreach my $l(@data) {
			chomp($l);
			next if($l=~/^\#/);
			my @t=split(/\s+/, $l);
			my $totalIntersect = scalar(computeClusterOverlap($t[14], $t[29]));
			print OUTFILE "$l\t$totalIntersect\n";
		}
		close OUTFILE;

		system("gzip $name{'rep1'}.bed");
		system("gzip $name{'rep2'}.bed");
		system("gzip $overlapFile");
	}
	## convert consitutively expressed loci file to BED, determine overlap between individual and all replicates
	elsif(defined($exprLociFile) && defined($replicates) && $replicates==2) {
		usage() if(!$overlapFile);

		## determine files corresponding to replicate 1 and both replicates
		my %file=(); my %name=();
		if($exprLociFile=~/Rep1/) {
			$file{'ind'}="$exprLociFile";
			$file{'all'}="$exprLociFile";
			$file{'all'}=~s/Rep1/All/;
			if(! -e $file{'ind'} || ! -e $file{'all'}) {
				print STDERR "Either of $file{'ind'} or $file{'all'} do not exist\n";
				exit(-1);
			}
			$name{'ind'}=$file{'ind'}; $name{'ind'}=~s/\.gz//g;
			$name{'all'}=$file{'all'}; $name{'all'}=~s/\.gz//g;
		}
		## determine files corresponding to replicate 2 and both replicates
		elsif($exprLociFile=~/Rep2/) {
			$file{'ind'}="$exprLociFile";
			$file{'all'}="$exprLociFile";
			$file{'all'}=~s/Rep2/All/;
			if(! -e $file{'ind'} || ! -e $file{'all'}) {
				print STDERR "Either of $file{'ind'} or $file{'all'} do not exist\n";
				exit(-1);
			}
			$name{'ind'}=$file{'ind'}; $name{'ind'}=~s/\.gz//g;
			$name{'all'}=$file{'all'}; $name{'all'}=~s/\.gz//g;
		}
		else { print STDERR "Cannot locate both replicates for $exprLociFile\n"; exit(-1); }

		## check format of individual replicate file
		my $found=`zless $file{'ind'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=(14+$tissues)) { print STDERR "Incorrect format for $file{'ind'}\n"; exit(-1); }

		## convert constitutively expressed loci file into bed format (replicate 1)
		system("echo \"#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\" > $name{'ind'}.bed");
		system("zgrep -v \"\^\\#\" $file{'ind'} | perl -an -F'/\\s+/' -e 'if(\$F[0]=~/bgP/) { \$strand=\"\+\"; } else { \$strand=\"\-\"; } \@t=split(/[\\:\\-]+/, \$F[1]); \$blocks=1; \$clusterInfo=(); \$uniqReads=100; foreach \$l(\@F[3..($tissues+2)]) { \@b=split(/\\|/,\$l); if(\$b[9]<\$uniqReads) { \$uniqReads=\$b[9]; } if(\$b[7]>\$blocks) { \$blocks=\$b[7]; } \$clusterInfo.=\"\$b[10],\" } \$clusterInfo=~s/\,\$//g; if(\$uniqReads>=$threshold{'minUniqReads'}) { print \"\$t[0]\\t\$t[1]\\t\$t[2]\\t\$F[0]\\t0\\t\$strand\\t\$F[$tissues+12]_\$F[$tissues+13]\\t\$F[$tissues+3]\\t\$F[$tissues+4]\\t\$F[$tissues+5]\\t\$F[$tissues+10]\\t\$F[2]\\t\$blocks\\t\$uniqReads\\t\$clusterInfo\\n\"; }' >> $name{'ind'}.bed");

		## check format of both replicate file
		$tissues=$tissues*2;
		$found=`zless $file{'all'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=(14+$tissues)) { print STDERR "Incorrect format for $file{'all'}\n"; exit(-1); }

		## convert constitutively expressed loci file into bed format (replicate 2)
		system("echo \"#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\" > $name{'all'}.bed");
		system("zgrep -v \"\^\\#\" $file{'all'} | perl -an -F'/\\s+/' -e 'if(\$F[0]=~/bgP/) { \$strand=\"\+\"; } else { \$strand=\"\-\"; } \@t=split(/[\\:\\-]+/, \$F[1]); \$blocks=1; \$clusterInfo=(); \$uniqReads=100; foreach \$l(\@F[3..($tissues+2)]) { \@b=split(/\\|/,\$l); if(\$b[9]<\$uniqReads) { \$uniqReads=\$b[9]; } if(\$b[7]>\$blocks) { \$blocks=\$b[7]; } \$clusterInfo.=\"\$b[10],\" } \$clusterInfo=~s/\,\$//g; if(\$uniqReads>=$threshold{'minUniqReads'}) { print \"\$t[0]\\t\$t[1]\\t\$t[2]\\t\$F[0]\\t0\\t\$strand\\t\$F[$tissues+12]_\$F[$tissues+13]\\t\$F[$tissues+3]\\t\$F[$tissues+4]\\t\$F[$tissues+5]\\t\$F[$tissues+10]\\t\$F[2]\\t\$blocks\\t\$uniqReads\\t\$clusterInfo\\n\"; }' >> $name{'all'}.bed");

		## determine overlapping coordinates between the two replicates
		my @data=();
		#@data = `bed2intersect.pl -f $file{'ind'}.bed -s $file{'all'}.bed | grep -v "^NA"`;
		@data = `intersectBed -a $name{'ind'}.bed -b $name{'all'}.bed -wo -f 0.80 -r -s | cut -f 1-30`;

		## determine overlap for 'pvclust' clusters between the two replicates
		open(OUTFILE,">$overlapFile");
		print OUTFILE "#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\tchr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\n";
		foreach my $l(@data) {
			chomp($l);
			next if($l=~/^\#/);
			print OUTFILE "$l\n";
		}
		close OUTFILE;

		system("gzip $name{'ind'}.bed");
		system("gzip $name{'all'}.bed");
		system("gzip $overlapFile");
	}
	## convert consitutively expressed loci file to BED
	elsif(defined($exprLociFile)) {
		## check format of second replicate file
		my $found=`zless $exprLociFile | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=(14+$tissues)) { print STDERR "Incorrect format for $exprLociFile\n"; exit(-1); }

		my $name=$exprLociFile; $name=~s/\.gz//g;

		## convert constitutively expressed loci file into bed format (replicate 1)
		system("echo \"#chr\tstart\tend\tid\tscore\tstrand\tanno\tsilhouette\tdunn\tclusterScore\tdbaScore\tclusters\tblocks\tuniqReads\tclusterInfo\" > $name.bed");
		system("zgrep -v \"\^\\#\" $exprLociFile | perl -an -F'/\\s+/' -e 'if(\$F[0]=~/bgP/) { \$strand=\"\+\"; } else { \$strand=\"\-\"; } \@t=split(/[\\:\\-]+/, \$F[1]); \$blocks=1; \$clusterInfo=(); \$uniqReads=100; foreach \$l(\@F[3..($tissues+2)]) { \@b=split(/\\|/,\$l); if(\$b[9]<\$uniqReads) { \$uniqReads=\$b[9]; } if(\$b[7]>\$blocks) { \$blocks=\$b[7]; } \$clusterInfo.=\"\$b[10],\" } \$clusterInfo=~s/\,\$//g; if(\$uniqReads>=$threshold{'minUniqReads'}) { if(!defined(\$uniqReads)) { \$uniqReads=100; } print \"\$t[0]\\t\$t[1]\\t\$t[2]\\t\$F[0]\\t0\\t\$strand\\t\$F[$tissues+12]_\$F[$tissues+13]\\t\$F[$tissues+3]\\t\$F[$tissues+4]\\t\$F[$tissues+5]\\t\$F[$tissues+10]\\t\$F[2]\\t\$blocks\\t\$uniqReads\\t\$clusterInfo\\n\"; }' >> $name.bed");
		system("gzip $name.bed");
	}
}
elsif($compute=~/[bB]+/) {
	## check, if window and step length is defined
	my $noSubWindow=0;
	if(!$win || !$steps) { $noSubWindow=1; }

	## open bed file
	my @data=();
	if(defined($bedFile)) { @data=openFile($bedFile); }
	else { @data=openFile(); }

	## header
	print "#Seq\tLength\tGC%\t";

	my $count=0;
	foreach my $l(@data) {
		if($l!~/^\#/ && $l!~/^$/) {
			my @BED=();
			parseBED($l, \@BED);

			## retrieve nucleotide sequence from UCSC and format it
			my $seq=`coor2seq.pl -c "$BED[0]:$BED[1]-$BED[2]" -s $BED[5] | grep -v "^>" | perl -ane '\$_=~tr/\\n/ /; \$_=~s/\\s+//g; print \$_;'`;
			## compute GC content of the sequence
			my $gc=0;
			for(split //, $seq) {
				$gc++ if($_=~/^G$/i || $_=~/^C$/i);
			}
			
			## compute nucleotide content for given coordinates
			my %frequency=();
			tie %frequency, 'Tie::IxHash';

			if($noSubWindow) { $win=length($seq); $steps=1; }

			## compute mononucleotide frequencies
			compNuclFreq($seq, $win, $steps, "mono", \%frequency);

			## compute dinucleotide frequencies
			compNuclFreq($seq, $win, $steps, "di", \%frequency);

			## compute trinucleotide frequencies
			compNuclFreq($seq, $win, $steps, "tri", \%frequency);

			## print the results
			if($count==0) {	foreach(keys(%frequency)) { print "$_\t"; } print "\n"; }
			print "$seq\t";
			print length($seq)."\t";
			printf("%0.2f\t", (($gc/length($seq))*100));
			foreach(keys(%frequency)) { print "$frequency{$_}\t"; } print "\n";

			$count++;
		}
	}
}
elsif($compute=~/[cC]+/) {
	usage() if(!$dpLociFile || !$deseqMatrixFile);

	## open candidate loci file
	my @data = openFile($dpLociFile);

	# header
	print "#coor\tfoldChange\tpval\tpadj\n";

	my $tissues=();
	foreach my $l(@data) {
		next if($l=~/^\#/);
		my @t=split(/\s+/, $l);

		if($l!~/^$/) {
			my %lociInfo=(); my %clusterInfo=();
			$lociInfo{'coor'} = "$t[0]:$t[1]-$t[2]";
			$lociInfo{'clusters'}=$t[11];

			my $i=1; my @clusters=split(/\,/,$t[14]);
			foreach(@clusters) {
				push(@{$clusterInfo{$_}}, $i);
				push(@{$lociInfo{'bgCluster'}}, $_);
				$i++;
			}

			my $deseqArg=();
			foreach my $cluster(sort { scalar(@{$clusterInfo{$a}}) <=> scalar(@{$clusterInfo{$b}}) } keys(%clusterInfo)) {
				if($lociInfo{'clusters'}>1 && $cluster!=0) {
					foreach(@{$lociInfo{'bgCluster'}}) {
						if($_==$cluster) { $deseqArg .= "DP,"; }
						else { $deseqArg .= "nDP,"; }
					} $deseqArg=~s/\,$//g;
					last;
				}
				elsif($lociInfo{'clusters'}==1) {
					foreach(@{$lociInfo{'bgCluster'}}) {
						if($_==$cluster) { $deseqArg .= "DP,"; }
						else { $deseqArg .= "nDP,"; }
					} $deseqArg=~s/\,$//g;
					last;
				}
			}
			#print "R --no-save --vanilla --slave < ~/software/myScripts/deseq.R --args $deseqMatrixFile 0 $deseqArg $lociInfo{'coor'} | grep -E \"^[0-9]+\"";
			my $deseq = `R --no-save --vanilla --slave < ~/software/myScripts/deseq.R --args $deseqMatrixFile 0 $deseqArg $lociInfo{'coor'} | grep -E "^[0-9]+"`;
			$deseq =~tr/\n/\t/;
			#print "$deseq\n";
			my @result = split(/\s+/, $deseq);
			print "$lociInfo{'coor'}\t$result[7]\t$result[8]\t$result[9]\n";
		}
	}
}
elsif($compute=~/^[dD]+$/) {
	my %bgInfo=(); my $first=(); my $second=();

	## open bed file
	my @data=();
	if(defined($bgFile)) { @data=openFile($bgFile); }
	else { @data=openFile(); }

	foreach(@data) {
		my @t=split(/\s+/, $_);
		if($_=~/^\>/) {
			if(%bgInfo) {
				#printf("%s\t%0.2f\t%0.2f\t%d\t%s\n", $bgInfo{'id'}, $first, $second, $bgInfo{'len'}, $bgInfo{'strand'});
				printf("%s\t%0.2f\t%0.2f\t%s\n", $bgInfo{'id'}, $first/($bgInfo{'len'}/2), $second/($bgInfo{'len'}/2), $bgInfo{'strand'}) if($bgInfo{'strand'}=~/\+/); 
				printf("%s\t%0.2f\t%0.2f\t%s\n", $bgInfo{'id'}, $second/($bgInfo{'len'}/2), $first/($bgInfo{'len'}/2), $bgInfo{'strand'}) if($bgInfo{'strand'}=~/\-/); 
			}
			$bgInfo{'id'}=$t[0];
			$bgInfo{'len'}=($t[3]-$t[2]);
			$bgInfo{'start'}=$t[2];
			$bgInfo{'mid'}=$t[2]+($bgInfo{'len'}/2);
			$bgInfo{'end'}=$t[3];
			$bgInfo{'strand'}=$t[4];
			$first=0; $second=0;
		}
		else {
			## compute read frequency in first half of the block group
			for(my $i=$bgInfo{'start'}; $i<$bgInfo{'mid'}; $i++) {
				if($i>=$t[1] && $i<=$t[2]) {
					$first+=$t[4];
				}
			}
			## compute read frequency in second half of the block group
			for(my $i=($bgInfo{'mid'}+1); $i<=$bgInfo{'end'}; $i++) {
				if($i>=$t[1] && $i<=$t[2]) {
					$second+=$t[4];
				}
			}
		}
	}
	if(%bgInfo) {
		#printf("%s\t%0.2f\t%0.2f\t%d\t%s\n", $bgInfo{'id'}, $first, $second, $bgInfo{'len'}, $bgInfo{'strand'});
		printf("%s\t%0.2f\t%0.2f\t%s\n", $bgInfo{'id'}, $first/($bgInfo{'len'}/2), $second/($bgInfo{'len'}/2), $bgInfo{'strand'}) if($bgInfo{'strand'}=~/\+/); 
		printf("%s\t%0.2f\t%0.2f\t%s\n", $bgInfo{'id'}, $second/($bgInfo{'len'}/2), $first/($bgInfo{'len'}/2), $bgInfo{'strand'}) if($bgInfo{'strand'}=~/\-/); 
	}
}
elsif($compute=~/^[eE]+$/) {
	$bins=10 if(!defined($bins));

	my %overlap=(); my $count=0; my $meanBinSize=();
	while(<STDIN>) {
		chomp($_);
		next if($_=~/^\#/ || $_=~/^$/);
		my @t=split(/\s+/, $_);
		my %bg=(); my %gene=(); my %bin=();
		tie %bin, 'Tie::IxHash';
		$bg{'start'}=$t[1];
		$bg{'end'}=$t[2];
		$bg{'loci'}=$t[3];
		$gene{'start'}=$t[7];
		$gene{'end'}=$t[8];
		$gene{'strand'}=$t[11];

		my $binSize=sprintf("%d", ($gene{'end'}-$gene{'start'})/$bins);

		my $j=1;
		$bin{1}{'start'}=$gene{'start'};
		while($j<$bins) {
			$bin{$j}{'end'}=$bin{$j}{'start'}+$binSize;
			$bin{$j}{'overlap'}=checkOverlap($bg{'start'}, $bg{'end'}, $bin{$j}{'start'}, $bin{$j}{'end'}, 3);
			#print STDERR "$bg{'start'}\t$bg{'end'}\t$bin{$j}{'start'}\t$bin{$j}{'end'}\t$bin{$j}{'overlap'}\n";
			$j++;
			$bin{$j}{'start'}=$bin{$j-1}{'end'};
		}
		$bin{$j}{'end'}=$gene{'end'};
		#print "$j\t".($gene{'end'}-$gene{'start'})."\n";
		$bin{$j}{'overlap'}=checkOverlap($bg{'start'}, $bg{'end'}, $bin{$j}{'start'}, $bin{$j}{'end'}, 3);
		#print STDERR "$bg{'start'}\t$bg{'end'}\t$bin{$j}{'start'}\t$bin{$j}{'end'}\t$bin{$j}{'overlap'}\n";

		my @keys=();
		if($gene{'strand'}=~/^\-$/) { @keys=reverse sort {$a <=> $b} keys(%bin); }
		else { @keys=sort {$a <=> $b} keys(%bin); }

		#print STDERR scalar(@keys)."\t$gene{'strand'}\n"; foreach(@keys) { print "$_\t"; } print "\n"; exit;
		for(my $i=0; $i<scalar(@keys); $i++) {
			$overlap{$i+1}+=$bin{$keys[$i]}{'overlap'};
			#if($bin{$keys[$i]}{'overlap'}>0) { $overlap{$i+1}++; }
			#print "$bin{$keys[$i]}{'overlap'}\t";
		}
		#print "$_\t$binSize\n";
		#print "\n";

		$meanBinSize+=$binSize;
		$count++;
	}

	foreach( sort {$a <=> $b } keys(%overlap)) {
		printf("%d\t%0.2f\t%0.2f\n", $_, $overlap{$_}/$count, $meanBinSize/$count);
	}
}
elsif($compute=~/[fF]+/) {
	usage() if(!$overlapFile || !$treeDir || !$sizeFactorFile);

	if(defined($id)) {
		if(! -d "$treeDir") { print STDERR "$treeDir does not exist\n"; exit; }
		my @data=openFile($overlapFile);

		my %dpInfo=();
		foreach my $l(@data) {
			next if($l=~/^\#/ || $l=~/^$/);
			next if($l!~/$id/);
			my @t=split(/\s+/,$l);
			$dpInfo{'rep1'}{'id'}=$id;
			$dpInfo{'rep2'}{'id'}="NA";
			$dpInfo{'rep1'}{'score'}=$t[$tissues+5];
			$dpInfo{'rep2'}{'score'}="NA";
			$dpInfo{'rep1'}{'clusterInfo'}="NA";
			$dpInfo{'rep2'}{'clusterInfo'}="NA";
			$dpInfo{'overlap'}="NA";

			my $pdfOut = "$dpInfo{'rep1'}{'id'}";
			next if (-f "$pdfOut.pdf");

			print STDERR "computing for $dpInfo{'rep1'}{'id'}\n";

			if($removeDummy) {
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep1 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep1 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep2 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep2");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep2 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep2");
			}
			else {
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep1 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep1 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep2 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep2");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep2 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep2");
			}
			print STDERR "blockbuster2ucsc (rep2).... done\n";

			$dpInfo{'rep1'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep1 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "multiTrack2pdf (rep1).... done\n";

			$dpInfo{'rep2'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep2 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep2");
			print STDERR "multiTrack2pdf (rep2).... done\n";

			#$dpInfo{'rep1'}{'pdfFile'}="bgN-1094802-1094847-chr10.conserved.norm.rep1.pdf";
			#$dpInfo{'rep2'}{'pdfFile'}="bgN-1094801-1094847-chr10.conserved.norm.rep2.pdf";
			chomp($dpInfo{'rep1'}{'pdfFile'}); chomp($dpInfo{'rep2'}{'pdfFile'});

			combinePDF($dpInfo{'rep1'}{'pdfFile'}, $dpInfo{'rep2'}{'pdfFile'}, $pdfOut, 2, %dpInfo);
			print STDERR "final output $pdfOut.. done\n";
		}
	}
	elsif(`head -n 1 $overlapFile | perl -ane 'print scalar(\@F)'`==32) {
		if(! -d "$treeDir/rep1") { print STDERR "$treeDir/rep1/ does not exist\n"; exit; }
		if(! -d "$treeDir/rep2") { print STDERR "$treeDir/rep2/ does not exist\n"; exit; }
		
		my @data=openFile($overlapFile);

		my %dpInfo=();
		foreach my $l(@data) {
			next if($l=~/^\#/ || $l=~/^$/);
			next if(defined($id) && $l!~/$id/);
			my @t=split(/\s+/,$l);
			#next if($t[11]==1 && $t[26]==1);
			#next if($t[30]<=5);
			$dpInfo{'rep1'}{'id'}=$t[3];
			$dpInfo{'rep2'}{'id'}=$t[18];
			$dpInfo{'rep1'}{'score'}=$t[9];
			$dpInfo{'rep2'}{'score'}=$t[24];
			$dpInfo{'rep1'}{'clusterInfo'}=$t[14];
			$dpInfo{'rep2'}{'clusterInfo'}=$t[29];
			$dpInfo{'overlap'}=$t[30];

			my $pdfOut = "$dpInfo{'rep1'}{'id'}_$dpInfo{'rep2'}{'id'}";
			next if (-f "$pdfOut.pdf");

			print STDERR "computing for $dpInfo{'rep1'}{'id'} and $dpInfo{'rep2'}{'id'}\n";

			system("blockbuster2ucsc.pl -b $treeDir/rep1/$dpInfo{'rep1'}{'id'}.norm -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
			system("blockbuster2ucsc.pl -b $treeDir/rep1/$dpInfo{'rep1'}{'id'}.norm -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

			system("blockbuster2ucsc.pl -b $treeDir/rep2/$dpInfo{'rep2'}{'id'}.norm -d -s > $dpInfo{'rep2'}{'id'}.norm.rep2");
			system("blockbuster2ucsc.pl -b $treeDir/rep2/$dpInfo{'rep2'}{'id'}.norm -w -s >> $dpInfo{'rep2'}{'id'}.norm.rep2");
			print STDERR "blockbuster2ucsc (rep2).... done\n";

			$dpInfo{'rep1'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep1 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "multiTrack2pdf (rep1).... done\n";

			$dpInfo{'rep2'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep2'}{'id'}.norm.rep2 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep2'}{'id'}.norm.rep2");
			print STDERR "multiTrack2pdf (rep2).... done\n";

			#$dpInfo{'rep1'}{'pdfFile'}="bgN-1094802-1094847-chr10.conserved.norm.rep1.pdf";
			#$dpInfo{'rep2'}{'pdfFile'}="bgN-1094801-1094847-chr10.conserved.norm.rep2.pdf";
			chomp($dpInfo{'rep1'}{'pdfFile'}); chomp($dpInfo{'rep2'}{'pdfFile'});

			combinePDF($dpInfo{'rep1'}{'pdfFile'}, $dpInfo{'rep2'}{'pdfFile'}, $pdfOut, 2, %dpInfo);
			print STDERR "final output $pdfOut.. done\n";
		}
	}
	elsif(`head -n 1 $overlapFile | perl -ane 'print scalar(\@F)'`==15) {
		if(! -d "$treeDir") { print STDERR "$treeDir does not exist\n"; exit; }
		my @data=openFile($overlapFile);

		my %dpInfo=();
		foreach my $l(@data) {
			next if($l=~/^\#/ || $l=~/^$/);
			next if(defined($id) && $l!~/$id/);
			my @t=split(/\s+/,$l);
			$dpInfo{'rep1'}{'id'}=$t[3];
			$dpInfo{'rep2'}{'id'}="NA";
			$dpInfo{'rep1'}{'score'}=$t[9];
			$dpInfo{'rep2'}{'score'}="NA";
			$dpInfo{'rep1'}{'clusterInfo'}=$t[14];
			$dpInfo{'rep2'}{'clusterInfo'}="NA";
			$dpInfo{'overlap'}="NA";

			my $pdfOut = "$dpInfo{'rep1'}{'id'}";
			next if (-f "$pdfOut.pdf");

			print STDERR "computing for $dpInfo{'rep1'}{'id'}\n";

			if($removeDummy) {
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep1 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep1 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep2 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep2");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | validateBlockbuster.pl -n 2 | grep Rep2 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep2");
			}
			else {
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep1 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep1 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep2 | blockbuster2ucsc.pl -d -s > $dpInfo{'rep1'}{'id'}.norm.rep2");
				system("subSampleBG.pl -f $treeDir/$dpInfo{'rep1'}{'id'}.norm -p | grep Rep2 | blockbuster2ucsc.pl -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep2");
			}
			print STDERR "blockbuster2ucsc (rep2).... done\n";

			$dpInfo{'rep1'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep1 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "multiTrack2pdf (rep1).... done\n";

			$dpInfo{'rep2'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep2 -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep2");
			print STDERR "multiTrack2pdf (rep2).... done\n";

			#$dpInfo{'rep1'}{'pdfFile'}="bgN-1094802-1094847-chr10.conserved.norm.rep1.pdf";
			#$dpInfo{'rep2'}{'pdfFile'}="bgN-1094801-1094847-chr10.conserved.norm.rep2.pdf";
			chomp($dpInfo{'rep1'}{'pdfFile'}); chomp($dpInfo{'rep2'}{'pdfFile'});

			combinePDF($dpInfo{'rep1'}{'pdfFile'}, $dpInfo{'rep2'}{'pdfFile'}, $pdfOut, 2, %dpInfo);
			print STDERR "final output $pdfOut.. done\n";
		}
	}
=cu
	elsif(`zgrep "^\\#" $overlapFile | perl -ane 'print scalar(\@F)'`==15) {
		if(! -d "$treeDir") { print STDERR "$treeDir does not exist\n"; exit; }
		my @data=openFile($overlapFile);

		my %dpInfo=();
		foreach my $l(@data) {
			next if($l=~/^\#/ || $l=~/^$/);
			my @t=split(/\s+/,$l);
			$dpInfo{'rep1'}{'id'}=$t[3];
			$dpInfo{'rep1'}{'score'}=$t[9];
			$dpInfo{'rep1'}{'clusterInfo'}=$t[14];

			print STDERR "computing for $dpInfo{'rep1'}{'id'}\n";

			system("blockbuster2ucsc.pl -b $treeDir/$dpInfo{'rep1'}{'id'}.norm -d -s > $dpInfo{'rep1'}{'id'}.norm.rep1");
			system("blockbuster2ucsc.pl -b $treeDir/$dpInfo{'rep1'}{'id'}.norm -w -s >> $dpInfo{'rep1'}{'id'}.norm.rep1");
			print STDERR "blockbuster2ucsc (rep1).... done\n";

			my $tracktitle="$dpInfo{'rep1'}{'id'}: $dpInfo{'rep1'}{'score'} ($dpInfo{'rep1'}{'clusterInfo'})";
			$dpInfo{'rep1'}{'pdfFile'}=`multiTrack2pdf.pl -t $dpInfo{'rep1'}{'id'}.norm.rep1 -l \"$tracktitle\" -s $sizeFactorFile | cut -f 2`;
			system("rm $dpInfo{'rep1'}{'id'}.norm.rep1");
			chomp($dpInfo{'rep1'}{'pdfFile'});
			print STDERR "multiTrack2pdf ($dpInfo{'rep1'}{'pdfFile'}).... done\n";
		}
	}
=cut
}
elsif($compute=~/[gG]+/) {
	usage() if(!$dpLociFile && !$overlapFile);

	## check format for differentially processed loci file
	if(defined($dpLociFile)) {
		my $found=`zless $dpLociFile | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=13) { print STDERR "Incorrect format for $dpLociFile\n"; exit(-1); }
	}

	if(defined($replicates) && $replicates==1) {
		## determine files corresponding to both replicate 1 and 2
		my %file=();
		if($dpLociFile=~/Rep1/) {
			$file{'rep1'}="$dpLociFile";
			$file{'rep2'}="$dpLociFile";
			$file{'rep2'}=~s/Rep1/Rep2/;
			if(! -e $file{'rep1'} || ! -e $file{'rep2'}) {
				print STDERR "Either of $file{'rep1'} or $file{'rep2'} do not exist\n";
				exit(-1);
			}
		}
		elsif($dpLociFile=~/Rep2/) {
			$file{'rep2'}="$dpLociFile";
			$file{'rep1'}="$dpLociFile";
			$file{'rep1'}=~s/Rep2/Rep1/;
			if(! -e $file{'rep1'} || ! -e $file{'rep2'}) {
				print STDERR "Either of $file{'rep1'} or $file{'rep2'} do not exist\n";
				exit(-1);
			}
		}
		else { print STDERR "Cannot locate both replicates for $dpLociFile\n"; exit(-1); }

		## check format of first replicate file
		my $found=`zless $file{'rep1'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=13) { print STDERR "Incorrect format for $file{'rep1'}\n"; exit(-1); }

		## check format of second replicate file
		$found=`zless $file{'rep2'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=13) { print STDERR "Incorrect format for $file{'rep2'}\n"; exit(-1); }

		my %sigLoci=();
		push(@{$sigLoci{'first'}},`zless $file{'rep1'} | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0]."\\n"; } }'`);
		push(@{$sigLoci{'second'}}, `zless $file{'rep2'} | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0]."\\n"; } }'`);

		## open overlap file
		my @data=openFile($overlapFile);

		## check, for differentially processed loci in overlap file
		my %dpStatus=();
		foreach my $l(@data) {
			chomp($l);
			if($l=~/^\#/) { print "$l\tdpStatusAcrossRep\n"; }
			my @t=split(/\s+/,$l);
			$dpStatus{'first'}=0;
			if(grep(/^$t[3]$/, @{$sigLoci{'first'}})) {
				$dpStatus{'first'}=1;
			}
			$dpStatus{'second'}=0;
			if(grep(/^$t[18]$/, @{$sigLoci{'second'}})) {
				$dpStatus{'second'}=1;
			}

			if($dpStatus{'first'}==1 && $dpStatus{'second'}==1) {
				print "$l\t0\n";
			}
			elsif($dpStatus{'first'}==1 && $dpStatus{'second'}==0) {
				print "$l\t1\n";
			}
			elsif($dpStatus{'first'}==0 && $dpStatus{'second'}==1) {
				print "$l\t2\n";
			}
		}
	}
	elsif(defined($replicates) && $replicates==2) {
		## determine files corresponding to replicate1 and all
		my %file=();
		if($dpLociFile=~/Rep1/) {
			$file{'ind'}="$dpLociFile";
			$file{'all'}="$dpLociFile";
			$file{'all'}=~s/Rep1/All/;
			if(! -e $file{'ind'} || ! -e $file{'all'}) {
				print STDERR "Either of $file{'ind'} or $file{'all'} do not exist\n";
				exit(-1);
			}
		}
		## determine files corresponding to replicate2 and all
		elsif($dpLociFile=~/Rep2/) {
			$file{'ind'}="$dpLociFile";
			$file{'all'}="$dpLociFile";
			$file{'all'}=~s/Rep2/All/;
			if(! -e $file{'ind'} || ! -e $file{'all'}) {
				print STDERR "Either of $file{'ind'} or $file{'all'} do not exist\n";
				exit(-1);
			}
		}
		else { print STDERR "Cannot locate both replicates for $dpLociFile\n"; exit(-1); }

		## check format of individual replicate file
		my $found=`zless $file{'ind'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=13) { print STDERR "Incorrect format for $file{'ind'}\n"; exit(-1); }

		## check format of all replicate file
		$found=`zless $file{'all'} | perl -ane 'print scalar(\@F); exit;'`;
		if($found!=13) { print STDERR "Incorrect format for $file{'all'}\n"; exit(-1); }

		my %sigLoci=();
		push(@{$sigLoci{'ind'}},`zless $file{'ind'} | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0]."\\n"; } }'`);
		$threshold{'minTissuesInCluster'}=4;
		push(@{$sigLoci{'all'}}, `zless $file{'all'} | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0]."\\n"; } }'`);

		## open overlap file
		my @data=openFile($overlapFile);

		## check, for differentially processed loci in overlap file
		my %dpStatus=();
		foreach my $l(@data) {
			chomp($l);
			if($l=~/^\#/) { print "$l\tdpStatusAcrossRep\n"; }
			my @t=split(/\s+/,$l);
			$dpStatus{'ind'}=0;
			if(grep(/^$t[3]$/, @{$sigLoci{'ind'}})) {
				$dpStatus{'ind'}=1;
			}
			$dpStatus{'all'}=0;
			if(grep(/^$t[18]$/, @{$sigLoci{'all'}})) {
				$dpStatus{'all'}=1;
			}

			if($dpStatus{'ind'}==1 && $dpStatus{'all'}==1) {
				print "$l\t0\n";
			}
			elsif($dpStatus{'ind'}==1 && $dpStatus{'all'}==0) {
				print "$l\t1\n";
			}
			elsif($dpStatus{'ind'}==0 && $dpStatus{'all'}==1) {
				print "$l\t2\n";
			}
		}
	}
	elsif(defined($replicates) && $replicates==3) {
		my %sigLoci=();

		push(@{$sigLoci{'all'}},`zless $dpLociFile | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { print \$F[0]."\\n"; }'`);

		## TO DEBUG / UNDERSTAND, USE THIS ONE LINER
		#zgrep -v "#" bgExpConsAll.stat | perl -ane '$sig=0; %clust=(); for($i=3; $i<=20; $i++) { $F[$i]=~s/\|.+\|/|/g; @t=split(/\|/,$F[$i]); $t[0]=~s/[0-9]+$//g; $clust{$t[1]}{$t[0]}++; } foreach $id(keys(%clust)) { print "$id\t"; $sum=0; foreach(keys(%{$clust{$id}})) { print "$_:$clust{$id}{$_}\t"; if($clust{$id}{$_}>=2 && $id!=0) { $sum+=$clust{$id}{$_}; } } if($sum>2) { $sig++; } print "\n"; } print "$sig\n\n";' | less

		## grep cluster information for overlapping loci from overlapFile
		system("zgrep \"^\#\" $overlapFile");
		foreach my $loci(@{$sigLoci{'all'}}) {
			chomp($loci);
			my $isSig=`zless $overlapFile | perl -ane 'if(\$F[0]=~/\^$loci\$/) { print \$_; }' | perl -ane '\$sig=0; %clust=(); for(\$i=3; \$i<=20; \$i++) { \$F[\$i]=~s/\\|.+\\|/|/g; \@t=split(/\\|/,\$F[\$i]); \$t[0]=~s/[0-9]+\$//g; \$clust{\$t[1]}{\$t[0]}++; } foreach \$id(keys(%clust)) { \$sum=0; foreach(keys(%{\$clust{\$id}})) { if(\$clust{\$id}{\$_}>=2 && \$id!=0) { \$sum+=\$clust{\$id}{\$_}; } } if(\$sum>2) { \$sig++; } } print \"\$sig\";'`;
			if($isSig=~/^$/) { print STDERR "Cannot find $loci in the first column of $overlapFile\n"; exit(-1); }
			if($isSig>=1) { system("zless $overlapFile | perl -ane 'if(\$F[0]=~/\^$loci\$/) { print \$_; }'"); }
		}
	}
	else {
		my %sigLoci=();
		#print "zless $dpLociFile | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0].\"\\n\"; } }'\n";
		push(@{$sigLoci{'all'}},`zless $dpLociFile | grep -v \"^\\#\" | perl -ane 'if(\$F[5]>=$threshold{'minClusterScore'} && \$F[9]>=$threshold{'minFScore'}) { %clust=(); foreach(\@{[split(/\\,/, \$F[12])]}) { \$clust{\$_}++ if(\$_>0); } \$max=0; foreach(keys(\%clust)) { if(\$clust{\$_}>\$max) { \$max=\$clust{\$_}; } } if(\$max>=$threshold{'minTissuesInCluster'}) { print \$F[0]."\\n"; } }'`);

		## grep cluster information for overlapping loci from overlapFile
		system("zgrep \"^\#\" $overlapFile");
		foreach my $loci(@{$sigLoci{'all'}}) {
			chomp($loci);
			#print("zless $overlapFile | perl -ane 'if(\$F[0]=~/\^$loci\$/) { print \$_; }'\n");
			system("zless $overlapFile | perl -ane 'if(\$F[0]=~/\^$loci\$/) { print \$_; }'");
		}
	}
}
elsif($compute=~/[hH]+/) {
	usage() if(!$treeDir);
	$overlapFileRep1All = "rep1_all.overlap" if(!$overlapFileRep1All);

	## retrieve file name for all replicate analysis corresponding to validated loci in replicate 1
	my %dpLoci=(); my @data=();
	if(defined($overlapFileRep1Rep2)) {	@data=openFile($overlapFileRep1Rep2);	}
	else { @data=openFile(); }
	foreach my $l(@data) {
		chomp($l);
		next if($l=~/^\#/);
		my @t=split(/\s+/,$l);
		$dpLoci{$t[3]}{'all_file'} = `zless $overlapFileRep1All | grep -w $t[3] | cut -f 19`;
		$dpLoci{$t[3]}{'des'} = "$l";
	}

	## determine deepBlockAlign scores for read profiles from same tissue across two replicates
	foreach my $loci(keys(%dpLoci)) {
		chomp($loci); chomp($dpLoci{$loci}{'all_file'});
		if(-e "$treeDir/$dpLoci{$loci}{'all_file'}.norm.aln") {
			$dpLoci{$loci}{'dba_score_same'} = `zless $treeDir/$dpLoci{$loci}{'all_file'}.norm.aln | perl -ane 'if(\$F[0]=~/Rep1/ && \$F[1]=~/Rep2/) { \$F[0]=~s/Rep1.*//g; \$F[1]=~s/Rep2.*//g; if(\$F[0]=~/^\$F[1]\$/) { \$t.="\$F[2],"; }} END { \$t=~s/\\,\$//; print "\$t"; }'`;
			$dpLoci{$loci}{'dba_score_diff'} = `zless $treeDir/$dpLoci{$loci}{'all_file'}.norm.aln | perl -ane '\$F[0]=~s/Rep.*//g; \$F[1]=~s/Rep.*//g; if(\$F[0]!~/^\$F[1]\$/) { \$t.="\$F[2],"; } END { \$t=~s/\\,\$//; print "\$t"; }'`;
		}
		else { print STDERR "Error: cannot find file, $treeDir/$dpLoci{$loci}{'all_file'}.norm.aln\n"; }
	}

	## print clustering and deepBlockAlign score information (same tissue)
	if(!defined($dbaScoreDiffTissue)) {
		foreach my $loci(keys(%dpLoci)) {
			my @t=split(/\s+/, $dpLoci{$loci}{'des'});
			my @sameClusterTissues = computeClusterOverlap($t[14], $t[29]);
			my @dbaScoreSame=split(/\,/,$dpLoci{$loci}{'dba_score_same'});

			#print "$dpLoci{$loci}{'des'}\n";
			for(my $i=0; $i<$tissues; $i++) {
				## same cluster
				if(grep(/$i/, @sameClusterTissues)) {
					print "$dbaScoreSame[$i]\t1\n";
				}
				## different cluster
				else {
					print "$dbaScoreSame[$i]\t0\n";
				}
			}	
			#print "$dpLoci{$loci}{'dba_score_same'}\n";
			#print "$dpLoci{$loci}{'dba_score_diff'}\n";
		}
	}
	## print deepBlockAlign score information (different tissue)
	else {
		foreach my $loci(keys(%dpLoci)) {
			my @dbaScoreDiff=split(/\,/,$dpLoci{$loci}{'dba_score_diff'});
			foreach(@dbaScoreDiff) {
				print "$_\n";
			}
		}
	}
}
elsif($compute=~/[iI]+/) {
	usage() if(!$bgExprLociFile);
	
	my @data=();
	if(defined($bgExprLociFile)) {
		@data=openFile($bgExprLociFile);
	}
	else { @data=openFile(); }

	my @file=();
	foreach my $l(@data) {
		my @t=split(/\s+/, $l);
		if($l=~/^\#/) { @file=@t;	}
		elsif($l!~/^$/) {
			my @coor=split(/\_/,$t[0]);
			my %clusterInfo=();
			for(my $i=1; $i<scalar(@t); $i++) {
				next if($t[$i]=~/^\-$/);
				if($t[$i]=~/\,/) {
					my @c=split(/\,/,$t[$i]);
					foreach(@c) {
						if($_=~/\(\+\)/) { $clusterInfo{'pos'}{'id'}=$_; $clusterInfo{'pos'}{'index'}=$i; }
						elsif($_=~/\(\-\)/) { $clusterInfo{'neg'}{'id'}=$_; $clusterInfo{'neg'}{'index'}=$i; }
					}	
				}
				else {
					if($t[$i]=~/\(\+\)/) { $clusterInfo{'pos'}{'id'}=$t[$i]; $clusterInfo{'pos'}{'index'}=$i; }
					elsif($t[$i]=~/\(\-\)/) { $clusterInfo{'neg'}{'id'}=$t[$i]; $clusterInfo{'neg'}{'index'}=$i; }
				}
			}
			if(defined($clusterInfo{'pos'})) {
				print "$coor[2]\t$coor[0]\t$coor[1]\tAll\t1\t+\t";
				$clusterInfo{'pos'}{'id'}=~s/\(\+\)//g;
                if(defined($bgDir)) {
    				my $loci=`zless $bgDir/$file[$clusterInfo{'pos'}{'index'}] | grep -m 1 -w $clusterInfo{'pos'}{'id'} | cut -f 11`;
	    			chomp($loci);
		    		print "$loci\n";
                }
                else { print "\n"; }
			}
			if(defined($clusterInfo{'neg'})) {
				print "$coor[2]\t$coor[0]\t$coor[1]\tAll\t1\t-\t";
				$clusterInfo{'neg'}{'id'}=~s/\(\-\)//g;
                if(defined($bgDir)) {
    				my $loci=`zless $bgDir/$file[$clusterInfo{'neg'}{'index'}] | grep -m 1 -w $clusterInfo{'neg'}{'id'} | cut -f 11`;
	    			chomp($loci);
		    		print "$loci\n";
                }
                else { print "\n"; }
			}
		}
	}
}
elsif($compute=~/[jJ]+/) {
	usage() if(!$dpLociFile || !$mapDir || !$treeDir);

	my @data=openFile($dpLociFile);

	foreach my $l(@data) {
		#next if($l!~/bgP_61641_61662_chr7.conserved/);
		next if($l=~/^\#/);
		my @t=split(/\s+/, $l);
		my %bgInfo=(); my @bgInfoNorm=();
		tie %bgInfo, 'Tie::IxHash';
		my $infile=();
		if($t[3]=~/bg[P|N]+\_/) {
			$infile="$t[3]";
		}
		elsif($t[0]=~/bg[P|N]+\_/) {
			$infile="$t[0]";
		}

		readbgInfo(\%bgInfo, $treeDir, "$infile.norm");

		## open output file
		chomp($outDir);
		next if(-e "$outDir/$infile.norm.readposbias");

		print STDERR "Normalizing read count for $infile.norm\n";
		open(OUTFILE, ">$outDir/$infile.norm.readposbias") || die $!;

		## check if the loci is composed of block groups having only one block
		my $is_one_block=0;
		if(`zgrep "^>" $treeDir/$infile | cut -f 8 | uniq | wc -l`==1) {
			$is_one_block=1;
		}

		foreach my $bg(keys(%bgInfo)) {
			push(@bgInfoNorm, normReadCountBySeqBias(\%{$bgInfo{$bg}}, $mapDir, $is_one_block));
		}

		## print normalized read count to the output file
		foreach(@bgInfoNorm) {
			print OUTFILE "$_\n";
		}

		close OUTFILE;
	}
}

##  subroutine to normalize read count with a block group using seqbias package (R)
sub normReadCountBySeqBias {
	my($bgInfo, $mapDir, $is_one_block) = @_;

	## compute bias
	my $yaml_file=$bgInfo->{'id'};
	$yaml_file=~s/\>//g;
	$yaml_file=~s/\_.+/.bam.yaml/g;

	#print "R --no-save --vanilla --slave < ~/software/myScripts/readPosBias.R --args 2 $mapDir/$yaml_file $refGenome \"$bgInfo->{'chr'}_$bgInfo->{'start'}_$bgInfo->{'end'}_$bgInfo->{'strand'}\"\n";
	my @bias=split(/\s+/, `R --no-save --vanilla --slave < ~/software/myScripts/readPosBias.R --args 2 $mapDir/$yaml_file $refGenome "$bgInfo->{'chr'}_$bgInfo->{'start'}_$bgInfo->{'end'}_$bgInfo->{'strand'}" 2>/dev/null`);
	#print "$yaml_file\t";
	#foreach(@bias) { print "$_\t"; } print "\n";
	#print scalar(@bias)."\t$bgInfo->{'start'}\t$bgInfo->{'end'}\n";

	## normalize read count
	my @temp=(); my $expr=();
	foreach my $b(sort { $a <=> $b } grep(/^[0-9]+$/, keys(%{$bgInfo}))) {
		foreach my $l(@{$bgInfo->{$b}{'info'}}) {
			my @t=split(/\s+/, $l);
			if($bgInfo->{'strand'}=~/^\+$/) {
				my $index=$t[1]-$bgInfo->{'start'};
				my $norm_count=sprintf("%0.2f", $t[4]/$bias[$index]);
				if($l!~/dummy/) {
					push(@temp, "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$norm_count\t$t[5]\t$t[6]");
					$expr+=$norm_count;
				}
				else {
					push(@temp, "$l");
					$expr+=$t[4];
				}
			}
			else {
				my $index=(scalar(@bias)-($bgInfo->{'end'}-$t[2]))-1;
				#print "$index\t".scalar(@bias)."\t$bgInfo->{'end'}\t$t[2]\n";
				my $norm_count=sprintf("%0.2f", $t[4]/$bias[$index]);
				if($l!~/dummy/) {
					push(@temp, "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$norm_count\t$t[5]\t$t[6]");
					$expr+=$norm_count;
				}
				else {
					push(@temp, "$l");
					$expr+=$t[4];
				}
			}
		}
	}

	## normalize read count of dummy block, if loci is composed of block groups having only one block
	if($is_one_block) {
		for(my $i=0; $i<scalar(@temp); $i++) {
			if($temp[$i]=~/dummy/) {
				my @t=split(/\s+/, $temp[$i]);
				$expr=$expr-$t[4];
				my $norm_count=sprintf("%0.2f", $expr/10);
				$temp[$i]="$t[0]\t$t[1]\t$t[2]\t$t[3]\t$norm_count\t$t[5]\t$t[6]";
				$expr=$expr+$norm_count;
			}
		}
	}

	my @bgInfoNorm=();
	push(@bgInfoNorm, "$bgInfo->{'id'}\t$bgInfo->{'chr'}\t$bgInfo->{'start'}\t$bgInfo->{'end'}\t$bgInfo->{'strand'}\t$expr\t$bgInfo->{'tags'}\t$bgInfo->{'blocks'}\t$bgInfo->{'name'}\t$bgInfo->{'anno'}\t$bgInfo->{'loci'}");
	push(@bgInfoNorm, @temp);
	return @bgInfoNorm;
}

## subroutine to reformat cluster overlap information
sub reformatClusterOverlap {
	my($clusters) = @_;

	my @t=split(/\,/, $clusters);

	my $i=0; my %seen=();
	foreach(@t) {
		if(!defined($seen{$_})) {
			$seen{$_}=$i;
			$i++;
		}
	}

	my $t=();
	foreach(@t) {
		$t.="$seen{$_},";
	}

	$t=~s/\,$//g;
	return $t;
}

## subroutine to compute cluster overlap between the two replicates
sub computeClusterOverlap {
	my ($clusterRep1, $clusterRep2) = @_;

	## organize cluster information
	my %clusterInfo=();

	## for replicate 1
	my $i=0;
	foreach(split(/\,/,$clusterRep1)) {
		push(@{$clusterInfo{'rep1'}{$_}}, $i);
		$i++;
	}
	#for replicate 2
	$i=0;
	foreach(split(/\,/,$clusterRep2)) {
		push(@{$clusterInfo{'rep2'}{$_}}, $i);
		$i++;
	}

	## compare clustering between the two replicates
	my @clusterSortRep1= reverse sort { scalar(@{$clusterInfo{'rep1'}{$a}}) <=> scalar(@{$clusterInfo{'rep1'}{$b}}) } keys(%{$clusterInfo{'rep1'}});
	my @clusterSortRep2=reverse sort { scalar(@{$clusterInfo{'rep2'}{$a}}) <=> scalar(@{$clusterInfo{'rep2'}{$b}}) } keys(%{$clusterInfo{'rep2'}});
	#print scalar(@clusterSortRep1)."\t".scalar(@clusterSortRep2)."\n";

	my $totalIntersect=0; my @sameClusterTissues=();
	for(my $i=0; $i<scalar(@clusterSortRep1); $i++) {
		my @maxIntersect=(); my $maxCluster=();
		#print scalar(@{$clusterInfo{'rep1'}{$clusterSortRep1[$i]}})."\t";
		for(my $j=0; $j<scalar(@clusterSortRep2); $j++) {
			#print "$clusterSortRep1[$i]\t$clusterSortRep2[$j]\n";
			my $lc=List::Compare->new(\@{$clusterInfo{'rep1'}{$clusterSortRep1[$i]}}, \@{$clusterInfo{'rep2'}{$clusterSortRep2[$j]}});
			my @intersection=$lc->get_intersection;
			#foreach(@{$clusterInfo{'rep1'}{$clusterSortRep1[$i]}}) { print "$_ "; } print "\n";
			#foreach(@{$clusterInfo{'rep2'}{$clusterSortRep2[$j]}}) { print "$_ "; } print "\n";
			#foreach(@intersection) { print "$_ "; } print "\n"; exit;
			if(scalar(@intersection)>scalar(@maxIntersect)) {
				@maxIntersect=@intersection;
				$maxCluster=$j;
			}
		}
		if(defined($maxCluster)) {
			#print scalar(@{$clusterInfo{'rep2'}{$clusterSortRep2[$maxCluster]}})."\t";
			#print scalar(@maxIntersect)."\t$maxCluster\n";
			$totalIntersect+=scalar(@maxIntersect);
			push(@sameClusterTissues, @maxIntersect);
			splice(@clusterSortRep2, $maxCluster, 1);
		}
	}

	#foreach(@sameClusterTissues) { print "$_ "; } print "\n";
	return @sameClusterTissues;
	#return $totalIntersect;
}
