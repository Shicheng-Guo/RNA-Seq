#!/usr/bin/perl -w

#use strict;
use warnings;
use Getopt::Long;
use perlModule;

##################################################################################################
use vars qw($fastaFile $minPercentIdentity $help);
$minPercentIdentity=0;

GetOptions ("f=s"  => \$fastaFile,
            "i=i"  => \$minPercentIdentity,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

##################################################################################################
sub usage {
  print STDERR "\nProgram: seqIdentify.pl (compute identities between all sequence pairs in the fasta file)\n";
  print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
  print STDERR "Version: 1.0\n";
  print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Credit: Stefan Seemann\n";
  print STDERR "Usage: seqIdentity.pl -f <file> [OPTIONS]\n";
  print STDERR " -f <file>       [fasta file | STDIN]\n";
  print STDERR "[OPTIONS]\n";
  print STDERR " -i <int>        [minimum percent identity (default: 0)]\n";
  print STDERR " -h <help>\n\n";
  exit(-1);
}

##################################################################################################

## open FASTA file
my @data=();
if(defined($fastaFile)) {
	@data=openFile($fastaFile);
}
else {
	@data=openFile();
}

$j=-1;
foreach(@data) {
	chomp $_;
	if(/^>/) {
		$j++;
		/>(.*)/;
		$name[$j] = $1;
	}
	else  {
		$_ =~ s/\./-/g;
		@tmp = split "";
		push @{$seq[$j]}, @tmp;
	}
}

$sum=0;
$nr=0;
for( $tj=0;$tj<=$j;$tj++ ) {
	for( $ttj=$tj+1;$ttj<=$j;$ttj++ ) {
		#compare $seq[$tj] with $seq[$ttj]
		$e=0;
		$l=$#{$seq[$tj]}+1;
		$updatedl = $l;
		for( $i=0;$i<$l;$i++ ) {
	    if( uc(${$seq[$tj]}[$i]) eq uc(${$seq[$ttj]}[$i]) ) {
				if( ${$seq[$tj]}[$i] eq "-" ) {
					$updatedl--;
				}
				else {
					$e++ 	
				}
    	}
		}
		if( $updatedl ) {
			$id = $e/$updatedl;
			print "$name[$tj]\t$name[$ttj]\t$id\t$e\t$updatedl\n" if($id>=$minPercentIdentity);
			$sum += $id;
			$nr++;
		}
	}
}

if($nr>0) { $id=$sum/$nr*100; }
else { $id=0; }
#print "average pairwise identity: $id\n";
print "$id\n";

