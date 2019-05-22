#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util;
use Cwd;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $segemehl $wigFile $log $title);
$title = "";

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions ("s=s"       => \$segemehl,
	    "o=s"       => \$wigFile,
	    "l"         => \$log,
	    "t=s"       => \$title,
	    "help"      => \$help,
             "h"        => \$help);
usage() if ($help || !$segemehl || !$wigFile);

my %wig = ();
readSegemehl($segemehl);
writeWIG("$wigFile", $title);



sub usage {
  print STDERR "\nusage: segemehl2wig.pl -s <file> -o <file>\n";
  print STDERR "write read density to a wig file\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -s <file>    segemehl output\n";
  print STDERR " -o <file>    wig files (file.pos.wig and file.neg.wig)\n";
  print STDERR "[OPTIONS]\n";
  print STDERR " -t <string>  title of the track ( default: (+) and (-) )\n";
  print STDERR " -l           logarithmize output ( default: off )\n";
  print STDERR " -h <file>    this (useful) help message\n";
  print STDERR "[VERSION]\n";
  print STDERR " 11-12-2010\n";
  print STDERR "[BUGS]\n";
  print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
  print STDERR "\n";
  exit(-1);
}


sub readSegemehl{
  my ($file) = @_;
  open(FILE, "<$file") || die "cannot open $file\n";
  while(<FILE>){
    chomp;
    next if(/\#/);
    my ($pairStat, $id, $dist, $qstart, $qend, $m, $mm, $in, $del, $strand, $start, $end, $chrom, $dum, $mapFreq) = split(/\s+/, $_);
    for(my $x = $start; $x <= $end; $x++){
      $wig{$strand}{$chrom}{$x} += 1 / $mapFreq;
    }
  }
  close(FILE);
}

sub writeWIG{
  my ($file, $title) = @_;

  open(POS, ">$file\.pos\.wig") || die "cannot open $file\.pos\.wig\n";
  print POS "track type=wiggle_0 name=\"$title (+)\" description=\"$title (+)\" visibility=full autoScale=off viewLimits=0.0:25.0 color=255,0,0\n";
  open(NEG, ">$file\.neg\.wig") || die "cannot open $file\.neg\.wig\n";
  print NEG "track type=wiggle_0 name=\"$title (-)\" description=\"$title (-)\" visibility=full autoScale=off viewLimits=0.0:25.0 color=0,0,255\n";

  foreach my $strand (keys %wig){
    foreach my $chrom (sort {$a cmp $b} keys %{$wig{$strand}}){
      if($strand eq "+") {print POS "variableStep chrom=$chrom\n";}
      if($strand eq "-") {print NEG "variableStep chrom=$chrom\n";}
      foreach my $p (sort {$a<=>$b} keys %{$wig{$strand}{$chrom}}){
	my $value = $wig{$strand}{$chrom}{$p};
	if($log){
	  $value = log($value + 1);
	}
	if($strand eq "+"){
          print POS "$p\t$value\n";
	}
	elsif($strand eq "-"){
	  print NEG "$p\t$value\n";
	}
      }
    }
  }
  close(POS);
  close(NEG);
}
