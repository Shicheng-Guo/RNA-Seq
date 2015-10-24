#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Tie::IxHash;
use perlModule;

############################################################################################################
# parse input options
use vars qw($fastq $adapter $length $help);

GetOptions ("i=s"  => \$fastq,
            "a=s"  => \$adapter,
            "l=i"  => \$length,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$adapter || !$length);

############################################################################################################
sub usage {
	print STDERR "\nProgram: clipAdapterWithBarCode.pl (clip adapter that have bar code from reads)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: clipAdapter.pl -i <file> -a <adapter> -l <int>\n";
	print STDERR " -i <file>       [first fasta file]\n";
	print STDERR " -j <file>       [second fasta file]\n";
	print STDERR " -a <adapter>    [adapter sequence]\n";
	print STDERR " -l <int>        [barcode length]\n\n";
	exit(-1);
}
############################################################################################################

my @data=openFile($fastq);

my %read=();
foreach(@data) {
	chomp($_);
	if($_=~/^\@/) {
		if(%read) {
			print "$read{'header'}\n";
			if(length($read{'seq'})<=18) {
				print "NNNNNNNNNNNNNNNNN\n";
				print "$read{'symbol'}\n";
				print "BBBBBBBBBBBBBBBBB\n";
			}
			else {
				print "$read{'seq'}\n";
				print "$read{'symbol'}\n";
				foreach(@{$read{'quality'}}[$read{'clipped'}..(scalar(@{$read{'quality'}})-1)]) {
					print "$_";
				}
				print "\n";
			}
		}
		%read=();
		$read{'header'}=$_;
	}
	elsif($_=~/^[ATGCN]+$/) {
		my $len1=length($_);

		## ADD ANY OTHER ADAPTER SEQUENCE HERE
		$_=~s/^.*AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT//g;

		$_=~s/^.*$adapter.{0,$length}//g;
		my $len2=length($_);
		$read{'seq'}=$_;
		$read{'clipped'}=$len1-$len2;
	}
	elsif($_=~/^\+/) {
		$read{'symbol'}=$_;
	}
	else {
		@{$read{'quality'}}=split(//,$_);
	}
}

if(%read) {
	print "$read{'header'}\n";
	print "$read{'seq'}\n";
	print "$read{'symbol'}\n";
	foreach(@{$read{'quality'}}[$read{'clipped'}..(scalar(@{$read{'quality'}})-1)]) {
		print "$_";
	}
	print "\n";
}
exit(0);
