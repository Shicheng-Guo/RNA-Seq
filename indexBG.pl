#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#########################################################################################
## parse input options
use vars qw($bgFile $help);

GetOptions ("i=s"  => \$bgFile,
						"help" => \$help,
            "h"    => \$help);

usage() if($help);

#########################################################################################
sub usage {
	print STDERR "\nProgram: indexBG.pl (create an index of block group file)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: indexBG.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>     [block group file | STDIN]\n";
	print STDERR " -h            [print this helpful message]\n\n";
	exit(-1);
}

#########################################################################################

## open text file
my $INFILE=();
if(defined($bgFile)) {
	if($bgFile=~/\.gz$/) { open($INFILE, "gunzip -c $bgFile |") || die $!; }
	else { open($INFILE, $bgFile) || die $!; }
}
else { $INFILE = *STDIN; }

my $i=0;
foreach my $l(<$INFILE>) {
	my @t=split(/\s+/, $l);
	if($l=~/^\>/) {
		$t[0]=~s/^\>//g;
		print "$t[0]\t$i\t".($i+$t[6])."\n";
	}
	$i++;
}

exit;
