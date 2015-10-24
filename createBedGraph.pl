#!/usr/bin/perl -w
##########################################################################################
## perl program to intersect two bed files
##########################################################################################

use strict;
use warnings;
use Getopt::Long;

#############################################################################################################
use vars qw($infile $bed $segemehl $blockbuster );
$overlapThreshold=1;

GetOptions ("f=s"  => \$bedFile1,
            "s=s"  => \$bedFile2,
            "o=i"  => \$overlapThreshold,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$bedFile1 || !$bedFile2);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: bed2anno.pl (retrieve annotations (bed) overlapping to query bed file)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: bed2intersect.pl -q <bed file> -a <bed file>\n";
	print STDERR " -f <query bed file>\n";
	print STDERR " -s <annotation bed file>\n";
	print STDERR " -o <overlap>          [default 1 nt]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}

#############################################################################################################
exit(0);
