#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

##########################################################################################################
# parse input options
my %coor=();
use vars qw($help) ;
$coor{'strand'}="any";

GetOptions ("c=s"  => \$coor{'coor'},
            "f=s"  => \$coor{'bgFile'},
            "s=s"  => \$coor{'strand'},
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$coor{'coor'} || !$coor{'bgFile'});

##########################################################################################################
sub usage {
	print STDERR "\nProgram: coor2bginfo.pl (given the coordinate, retrieve block group information)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2bginfo.pl -c <coordinate> -f <file> [OPTIONS]\n";
	print STDERR " -c <coordinate>  [chr:start-end]\n";
	print STDERR " -f <file>        [block group file in blockbuster format]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -s <strand>      [+/- (default: any)]\n\n";
	exit(-1);
}
##########################################################################################################

## initialize strand
if($coor{'strand'}=~/any/) { $coor{'strand'}=".+"; }
else { $coor{'strand'} = "\\$coor{'strand'}"; }

($coor{'chr'}, $coor{'start'}, $coor{'end'}) = split(/[\:\-]+/, $coor{'coor'});

if(!$coor{'chr'} || !$coor{'start'} || !$coor{'end'}) { print STDERR "Incorrect coordinate\n"; usage(); }

$coor{'start'}=~s/\,//g;
$coor{'end'}=~s/\,//g;

my @data=();
if($coor{'bgFile'}=~/\.gz$/) {
	@data=`zcat $coor{'bgFile'}`;
}
else {
	open(INFILE, $coor{'bgFile'}) || die $!;
	@data=<INFILE>;
	close INFILE;
}

foreach my $l(@data) {
	next if($l!~/^\>/);
	my @bgInfo = split(/\s+/, $l);
	next if($bgInfo[1]!~/^$coor{'chr'}$/);
	next if($bgInfo[4]!~/^$coor{'strand'}$/);

	if($bgInfo[2]>=$coor{'start'} && $bgInfo[3]<=$coor{'end'} ||
	$bgInfo[2]>=$coor{'start'} && $bgInfo[2]<=$coor{'end'} ||
	$bgInfo[3]>=$coor{'start'} && $bgInfo[3]<=$coor{'end'} ||
	$bgInfo[2]<=$coor{'start'} && $bgInfo[3]>=$coor{'end'}) {
		#print "$bgInfo[0]\t$bgInfo[1]\t$bgInfo[2]\t$bgInfo[3]_$bgInfo[6]\t$bgInfo[4]\t$bgInfo[5]\n";
		$bgInfo[0]=~s/\>//;
		#print "$coor{'coor'}\t$l\n";
		my $bgCount=`zless $coor{'bgFile'} | egrep "\\<$bgInfo[0]\\>.*\\<$bgInfo[2]\\>.*\\<$bgInfo[3]\\>" | wc -l`;
		if($bgCount==1) {
				system("zless $coor{'bgFile'} | egrep \"\\<$bgInfo[0]\\>.*\\<$bgInfo[2]\\>.*\\<$bgInfo[3]\\>\" -A $bgInfo[6]");
			}
		else {
			print STDERR "Error: more than one occurence of $bgInfo[0] in $coor{'bgFile'}\n";
			print STDERR "Error: zless $coor{'bgFile'} | egrep \"\\\<$bgInfo[0]\\\>\.*\\\<$bgInfo[2]\\\>\.*\\\<$bgInfo[3]\\\>\"\n";
			exit(-1);
		}
	}
}

exit;
