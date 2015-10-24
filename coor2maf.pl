#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

##################################################################################################
use vars qw($bedFile $maf_location $outDir $help);
$outDir=".";
$maf_location="/chromo/ucsc/gbdb/hg19/multiz46way/maf";

GetOptions ("i=s"  => \$bedFile,
            "m=s"  => \$maf_location,
            "o=s"  => \$outDir,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$outDir);

##################################################################################################
sub usage {
	print STDERR "\nProgram: coor2maf.pl (retrieve maf blocks corresponding to the coordinates in BED format)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: coor2maf.pl -i <file> -m <dir>\n";
	print STDERR " -i <file>        [coordinates in BED format | STDIN]\n";
	print STDERR " -m <dir>         [directory or file having MAF blocks (default: /chromo/ucsc/gbdb/hg19/multiz46way/maf)]\n";
	print STDERR " -o <dir>         [output directory (default: current)]\n";
	print STDERR " -h <help>\n\n";
	exit(-1);
}

##################################################################################################

## open coordinate file
my @data=();
if(defined($bedFile)) {
	@data=openFile($bedFile);
}
else {
	@data=openFile();
}

foreach my $l(@data) {
	chomp($l);
	my($chr,$start,$end,$name,$score,$strand) = split(/\s+/, $l);
	if($start!~/^[0-9]+$/ || $end!~/^[0-9]+$/) {
		print STDERR "Incorrect BED format\n\t$l\n"; exit;
	}
	if(defined($name)) { $name=~s/[^A-Za-z0-9]+//g; }
	else { $name="$chr"."_"."$start"."_"."$end"; }
	my $tmpName="$chr"."_"."$start"."_"."$end"."_"."$name";
	system("echo \'$l\' > $outDir/$tmpName.tmp");

	## retrieve MAF block corresponding to the coordinate
	if(-d $maf_location) {
		opendir(INDIR, $maf_location) || die $!;
		my @file_list = grep(/[\.\_]*$chr[\.\_]*/, readdir(INDIR));

		my $maf_files=();
		foreach(@file_list) { $maf_files.="$maf_location/$_ "; }
		#print "$tmpName\t$maf_files\n"; exit;
		system("mafsInRegion $outDir/$tmpName.tmp $outDir/$tmpName.maf $maf_files > /dev/null");
		close INDIR;
	}	
	else {
		system("mafsInRegion $outDir/$tmpName.tmp $outDir/$tmpName.maf $maf_location > /dev/null");
	}
	system("rm $outDir/$tmpName.tmp");
}

exit(0);
