#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $blockbuster $bedFile $wigFile $log $title $sepTitle $stdout);
$title = "";
$sepTitle=0;

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("b=s"  => \$blockbuster,
            "o=s"  => \$wigFile,
            "l"    => \$log,
            "t=s"  => \$title,
            "s=s"  => \$sepTitle,
            "u"    => \$stdout,
            "help" => \$help,
            "h"    => \$help);
usage() if ($help || !$blockbuster || (!$wigFile && !$stdout));

my (%pos, %neg, %blocks) = ();

open(BEDFILE, ">$wigFile\.bed") || die "cannot open $wigFile\.bed\n" if(!$stdout);
open(POSWIGFILE, ">$wigFile\.pos\.wig") || die "cannot open $wigFile\.pos\.wig\n" if(!$stdout);
open(NEGWIGFILE, ">$wigFile\.neg\.wig") || die "cannot open $wigFile\.neg\.wig\n" if(!$stdout);
if($sepTitle==0) {
	print "track name=\"$title - blockbuster blocks\" description=\"$title - blockbuster blocks\"  itemRgb=\"On\" visibility=\"3\"\n" if($stdout);
	print BEDFILE "track name=\"$title - blockbuster blocks\" description=\"$title - blockbuster blocks\"  itemRgb=\"On\" visibility=\"3\"\n" if(!$stdout);
	print POSWIGFILE "track type=wiggle_0 name=\"$title (+)\" description=\"$title (+)\" visibility=full autoScale=on viewLimits=0.0:25.0 color=255,0,0\n" if(!$stdout);
	print NEGWIGFILE "track type=wiggle_0 name=\"$title (-)\" description=\"$title (-)\" visibility=full autoScale=on viewLimits=0.0:25.0 color=0,0,255\n" if(!$stdout);
}

readClusters($blockbuster);

close BEDFILE if(!$stdout);
close POSWIGFILE if(!$stdout);
close NEGWIGFILE if(!$stdout);

sub usage {
	print STDERR "\nusage: blockbuster2ucsc.pl -b <file> -o <file>\n";
	print STDERR "write read density to a wig file and blocks to a bed file\n";
	print STDERR "\n";
	print STDERR "[INPUT]\n";
	print STDERR " -b <file>    blockbuster output\n";
	print STDERR " -o <file>    output files (file.bed, file.pos.wig and file.neg.wig)\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -t <string>  title of the track ( default: (+) and (-) )\n";
	print STDERR " -l           logarithmize output ( default: off )\n";
	print STDERR " -s <int>     seperate track title for each block group ( default: (0) and (1) )\n";
	print STDERR " -u           output to STDOUT in bed format\n";
	print STDERR " -h <file>    this (useful) help message\n";
	print STDERR "[VERSION]\n";
	print STDERR " 10-10-2012\n";
	print STDERR "[BUGS]\n";
	print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
	print STDERR "\n";
	exit(-1);
}


sub readClusters{
	my ($file) = @_;
	open(INFILE, "<$file") || die "cannot open $file\n";
	my ($clusterID, $chrom, $clusterStart, $clusterEnd, $strand, $clusterExp, $clusterReadCound, $clusterBlockCount, $ncRNAID, $ncRNAtype, $ncRNAclass, $dum, $start, $end, $id, $exp, $block) = ();
	while(<INFILE>) {
		chomp;
		next if(/\#/);
		if(/^\>/) {
			($clusterID, $chrom, $clusterStart, $clusterEnd, $strand, $clusterExp, $clusterReadCound, $clusterBlockCount, $ncRNAID, $ncRNAtype, $ncRNAclass) = split(/\s+/, $_);
			if(%blocks) {
				writeBED($title);
				writeWIG($title) if(!$stdout);
				(%pos, %neg, %blocks) = ();
			}
			$title=$clusterID if($sepTitle);
		}
		else {
			($chrom, $start, $end, $id, $exp, $strand, $block) = split(/\s+/, $_);
			if(!exists($blocks{$clusterID}{$block})) {
				$blocks{$clusterID}{$block}{start} = $start;
				$blocks{$clusterID}{$block}{end} = $end;
				$blocks{$clusterID}{$block}{strand} = $strand;
				$blocks{$clusterID}{$block}{chrom} = $chrom;
				$blocks{$clusterID}{$block}{expr} = $exp;
			} else {
				$blocks{$clusterID}{$block}{expr} += $exp;
				if($start < $blocks{$clusterID}{$block}{start}) {
					$blocks{$clusterID}{$block}{start} = $start;
				}
				if($end > $blocks{$clusterID}{$block}{end}){
					$blocks{$clusterID}{$block}{end} = $end;
				}
			}
			for(my $x = $start; $x <= $end; $x++) {
				if($strand eq "+") {
					$pos{$chrom}{$x} += $exp;
				}
				else {
					$neg{$chrom}{$x} += $exp;
				}
			}
		}
	}
	if(%blocks) {
		$title=$clusterID if($sepTitle);
		writeBED($title);
		writeWIG($title) if(!$stdout);
		(%pos, %neg, %blocks) = ();
	}
	close(INFILE);
}

sub writeBED{
	my ($title) = @_;
	print  "track name=\"$title - blockbuster blocks\" description=\"$title - blockbuster blocks\"  itemRgb=\"On\" visibility=\"3\"\n" if($sepTitle && $stdout);
	print  BEDFILE "track name=\"$title - blockbuster blocks\" description=\"$title - blockbuster blocks\"  itemRgb=\"On\" visibility=\"3\"\n" if($sepTitle && !$stdout);
	foreach my $cluster (keys %blocks) {
		foreach my $block (keys %{$blocks{$cluster}}) {
			print "$blocks{$cluster}{$block}{chrom}\t$blocks{$cluster}{$block}{start}\t$blocks{$cluster}{$block}{end}\t$cluster"."_"."$blocks{$cluster}{$block}{expr}\t0\t$blocks{$cluster}{$block}{strand}\t$blocks{$cluster}{$block}{start}\t$blocks{$cluster}{$block}{end}\t255,0,0\n" if($stdout);
			print BEDFILE "$blocks{$cluster}{$block}{chrom}\t$blocks{$cluster}{$block}{start}\t$blocks{$cluster}{$block}{end}\t$cluster"."_"."$blocks{$cluster}{$block}{expr}\t0\t$blocks{$cluster}{$block}{strand}\t$blocks{$cluster}{$block}{start}\t$blocks{$cluster}{$block}{end}\t255,0,0\n" if(!$stdout);
		}
	}
}

sub writeWIG{
	my ($title) = @_;

	foreach my $chrom (sort {$a cmp $b} keys %pos) {
		print POSWIGFILE "track type=wiggle_0 name=\"$title (+)\" description=\"$title (+)\" visibility=full autoScale=on viewLimits=0.0:25.0 color=255,0,0\n" if($sepTitle);
		print POSWIGFILE "variableStep chrom=$chrom\n";
		foreach my $p (sort {$a<=>$b} keys %{$pos{$chrom}}) {
			my $value = $pos{$chrom}{$p};
			if($log){ $value = log($value + 1); }
			print POSWIGFILE "$p\t$value\n";
		}
	}

	foreach my $chrom (sort {$a cmp $b} keys %neg) {
		print NEGWIGFILE "track type=wiggle_0 name=\"$title (-)\" description=\"$title (-)\" visibility=full autoScale=on viewLimits=0.0:25.0 color=0,0,255\n" if($sepTitle);
		print NEGWIGFILE "variableStep chrom=$chrom\n";
		foreach my $p (sort {$a<=>$b} keys %{$neg{$chrom}}) {
			my $value = $neg{$chrom}{$p};
			if($log){ $value = log($value + 1); }
			print NEGWIGFILE "$p\t$value\n";
		}
	}
}
