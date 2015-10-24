#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

##########################################################################################################
use vars qw($ucscTrackFile $colString $help);

my %color=(0, "216,222,174",
           1, "47,47,47",
           2, "204,182,71",
           3, "206,120,152",
           4, "152,197,171",
           5, "39,66,87",
           6, "214,222,228");

GetOptions ("u=s"  => \$ucscTrackFile,
            "c=s"  => \$colString,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$colString);

##########################################################################################################
sub usage {
	print STDERR "\nProgram: colorUCSCTrack.pl (color UCSC tracks)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: colorUCSCTrack.pl -u <file> -c <string>\n";
	print STDERR " -u <file>    [file with UCSC tracks | STDIN]\n";
	print STDERR " -c <string>  [string with color information (series of int seperated by comma)]\n";
	print STDERR " -h <help>    [this useful help message]\n\n";
	exit(-1);
}

#########################################################################################################

## open UCSC track file
my $INFILE=(); my @data=();
if(defined($ucscTrackFile)) {
	if($ucscTrackFile=~/\.gz$/) { open($INFILE, "gunzip -c $ucscTrackFile |") || die $!; }
	else { open($INFILE, $ucscTrackFile) || die $!; }
	@data=<$INFILE>;
	close $INFILE;
}
else { $INFILE=*STDIN; @data=<$INFILE>; }

## parse the color string
my @colorAttribute=split(/\,/, $colString);

my $TRACK_COUNTER=0;
foreach my $l(@data) {
	chomp($l);
	if($l=~/track type\=bedGraph/) {
		if(defined($colorAttribute[$TRACK_COUNTER])) { 
			print "$l color=$color{$colorAttribute[$TRACK_COUNTER]}\n";
		}
		else {
			print "$l color=$color{6}\n";
		}
		$TRACK_COUNTER++;
	}
	else { print "$l\n"; }
}

exit;
