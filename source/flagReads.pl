#!/usr/bin/perl -w
#usage: perl flagReads.pl <anno file> <map file> 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use vars qw ($help $annoFile $mapFile);

GetOptions ("a=s"  => \$annoFile,
	    "m=s"  => \$mapFile,
	    "help" => \$help,
	    "h"    => \$help);

usage() if($help || !$annoFile || !$mapFile);

sub usage {
	print STDERR "\nusage: flagReads.pl -a <file> -m <file> [OPTIONS]\n";
	print STDERR "flag reads overlapping with annotated ncRNAs or loci\n";
	print STDERR "\n";
	print STDERR "[INPUT]\n";
	print STDERR " -c <file>      annotation file\n";
	print STDERR " -a <file>      segemehl map file\n";
	print STDERR "[VERSION]\n";
	print STDERR " 11-29-2011\n";
	print STDERR " Please report bugs to sachin\@rth.dk\n";
	print STDERR "\n";
	exit(-1);
}

open(INFILE, $annoFile) || die $!;

my %anno=();
foreach(<INFILE>) {
	my @t = split(/\s+/, $_);
	$anno{"$t[0]_$t[1]_$t[2]"}{'chr'} = $t[0];
	$anno{"$t[0]_$t[1]_$t[2]"}{'start'} = $t[1];
	$anno{"$t[0]_$t[1]_$t[2]"}{'end'} = $t[2];
	$anno{"$t[0]_$t[1]_$t[2]"}{'anno'} = $t[3];
	$anno{"$t[0]_$t[1]_$t[2]"}{'strand'} = $t[5];
}

close INFILE;

#%anno=(); $anno{"chr1_12590391_12590475"}{'chr'} = "chr1"; $anno{"chr1_12590391_12590475"}{'start'} = 12590391; $anno{"chr1_12590391_12590475"}{'end'} = 12590475; $anno{"chr1_12590391_12590475"}{'anno'} = "hY3_RNA-related"; $anno{"chr1_12590391_12590475"}{'strand'} = "-";

my @mapData=(); my $found=0; my %seen=(); my @reads=();

foreach my $key(sort { $anno{$a}{'chr'} cmp $anno{$b}{'chr'} } keys(%anno)) {
	print ">$key\t$anno{$key}{'chr'}\t$anno{$key}{'start'}\t$anno{$key}{'end'}\t$anno{$key}{'anno'}\t$anno{$key}{'strand'}";
	$found=0; @reads=();
	if(!defined($seen{$anno{$key}{'chr'}})) { @mapData = `zcat $mapFile | grep -w $anno{$key}{'chr'}`; $seen{$anno{$key}{'chr'}}=1; }

	foreach(@mapData) {
		chomp($_);
		my @t = split(/\s+/, $_);
		if($t[10]>=$anno{$key}{'start'} && $t[11]<=$anno{$key}{'end'} && $t[12]=~/^$anno{$key}{'chr'}$/ && $t[9]=~/^\Q$anno{$key}{'strand'}\E$/) {
			push(@reads, $_);
			$found=1;
		}
		elsif($t[10]>=$anno{$key}{'start'} && $t[11]<=$anno{$key}{'end'} && $t[12]=~/^$anno{$key}{'chr'}$/) {
			push(@reads, $_);
			$found=2;
		}
	}
	if($found==0) { print "\t0\t0\n"; }
	elsif($found==1) { print "\t1\t".scalar(@reads)."\n"; foreach(@reads) { print "$_\n"; } }
	elsif($found==2) { print "\t2\t".scalar(@reads)."\n"; foreach(@reads) { print "$_\n"; } }
}

exit;
