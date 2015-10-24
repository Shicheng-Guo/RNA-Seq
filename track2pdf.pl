#!/usr/bin/perl -w

use strict;
use warnings;

use GD::Simple;
use Getopt::Long;
use POSIX qw(strftime);

############################################################################################################
## parse input options
use vars qw($bedfile $trackfile $tmpFilePrefix $sizeFactor $png $help);
$tmpFilePrefix="";

GetOptions ("b=s"  => \$bedfile,
	          "t=s"  => \$trackfile,
            "x=s"  => \$tmpFilePrefix,
            "s=s"  => \$sizeFactor,
            "p"    => \$png,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$trackfile);

############################################################################################################
sub usage {
	print STDERR "\nProgram: track2pdf.pl (visualize ucsc track file correponding to BED coordinates)\n";
	print STDERR "Author: University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: track2pdf.pl -b <file> -t <file> -p <file> [OPTIONS]\n";
	print STDERR " -b <file>     [input BED file | STDIN]\n";
	print STDERR " -t <file>     [input UCSC track file]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -x <string>   [temporary prefix for output pdf file]\n";
	print STDERR " -s <float>    [size factor to normalize the block expression]\n";
	print STDERR " -p            [show track as png image (default: pdf)]\n";
	print STDERR " -h            [help message]\n";
	exit(-1);
}
############################################################################################################

my %data=();
## read BED file
if(defined($bedfile)) {
	if($bedfile=~/\.gz/) { open(INFILE, "gunzip -c $bedfile |") || die $!; }
	else { open(INFILE, $bedfile) || die $!; }
	@{$data{'bed'}}=<INFILE>;
	close INFILE;
}
else {
	my $INFILE = *STDIN;
	@{$data{'bed'}}=<$INFILE>;
}

## read UCSC track file
if($trackfile=~/\.gz/) { open(INFILE, "gunzip -c $trackfile |") || die $!; }
else { open(INFILE, $trackfile) || die $!; }
@{$data{'track'}}=<INFILE>;
close INFILE;

my $count=1;
foreach my $d(@{$data{'bed'}}) {
	chomp($d);
	my @bed_coor=split(/\s+/, $d);
	if(scalar(@bed_coor)<5) {
		print STDERR "Error: incorrect BED format at \"$d\"\n";
		exit(-1);
	}
	my %image=('start'=>1000000000000000000000000, 'end'=>0, 'maxExpr'=>0);
	## determine chromosome
	$image{'chr'}=$bed_coor[0];
	## determine strand orientation
	$image{'strand'}=$bed_coor[5];

	my %track=(); my $track_name=(); my $track_strand=(); my $track_chr=();
	foreach my $t(@{$data{'track'}}) {
		my @track_coor=split(/\s+/, $t);
		if($t=~/name\=.*blockbuster blocks/) {
			$track_name="blockbuster";
			if($t=~/\(\+\)/) { $track_strand="+"; }
			elsif($t=~/\(\-\)/) { $track_strand="-"; }
		}
		elsif($t=~/name\=.*bigWig/) {
			$track_name="bigWig";
			if($t=~/\(\+\)/) { $track_strand="+"; }
			elsif($t=~/\(\-\)/) { $track_strand="-"; }
		}
		elsif($t=~/name\=.*wig/) {
			$track_name="wig";
			if($t=~/\(\+\)/) { $track_strand="+"; }
			elsif($t=~/\(\-\)/) { $track_strand="-"; }
		}
		elsif($t=~/variableStep/) {
			$track_chr=$t;
			$track_chr=~s/^.+\=//g;
		}
		elsif(defined($track_name) && $t!~/^$/ && $track_coor[0]=~/^$bed_coor[0]$/ && $track_coor[1]>=$bed_coor[1] && $track_coor[2]<=$bed_coor[2]) {
			next if(defined($track_strand) && $track_strand!~/\Q$image{'strand'}\E/);
			if($track_name=~/blockbuster/) { 
				push(@{$track{$track_name}}, $t);
			}
			elsif($track_name=~/bigWig/) {
				push(@{$track{$track_name}}, $t);
			}
		}
		elsif(defined($track_name) && $track_name=~/wig/ && $track_chr=~/^$bed_coor[0]$/ && $track_coor[0]>=$bed_coor[1] && $track_coor[0]<=$bed_coor[2]) {
			next if(defined($track_strand) && $track_strand!~/\Q$image{'strand'}\E/);
			push(@{$track{$track_name}}, $t);
		}
	}

	## nullify wig track, if both bigWig and wig are present
	@{$track{'wig'}}=() if(defined($track{'wig'}) && defined($track{'bigWig'}));

	## error, if all tracks are empty
	if(!defined($track{'blockbuster'}) && !defined($track{'bigWig'}) && !defined($track{'wig'})) {
		print STDERR "Error: no track observed for \"$d\" in $trackfile\n";
		exit(-1);
	}

	## determine range of image
	foreach(@{$track{'blockbuster'}}) {
		my @t=split(/\s+/,$_);
		if($t[1]<$image{'start'}) { $image{'start'}=$t[1]; }
		if($t[2]>$image{'end'}) { $image{'end'}=$t[2]; }
	}
	## determine maximum expression (bigWig)
	foreach(@{$track{'bigWig'}}) {
		my @t=split(/\s+/,$_);
		if($t[3]>$image{'maxExpr'}) { $image{'maxExpr'}=$t[3]; }
	}
	## determine maximum expression (wig)
	foreach(@{$track{'wig'}}) {
		my @t=split(/\s+/,$_);
		if($t[1]>$image{'maxExpr'}) { $image{'maxExpr'}=$t[1]; }
	}

	## image width
	$image{'width'}=($image{'end'}-$image{'start'})+1;

	## open pdf file for writing
	my $pdffile="$bed_coor[0]_$bed_coor[1]_$bed_coor[2]_$bed_coor[3]_$bed_coor[5]";
	open(OUTFILE, ">$tmpFilePrefix"."$pdffile.tex") || die $!;

print OUTFILE <<LATEX;
\\documentclass{standalone}
\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{tikz}
\\usepackage{subfigure}
\\usepackage{adjustbox}
\\definecolor{dark-red}{RGB}{100,0,0}
\\definecolor{light-red}{RGB}{70,0,0}
\\definecolor{dark-green}{RGB}{0,100,0}
\\definecolor{dark-blue}{RGB}{0,0,100}
\\definecolor{black}{RGB}{0,0,0}
\\begin{document}
\\centering
\\begin{adjustbox}{width=15 cm,height=5cm}
\\begin{tikzpicture}[scale=1]
LATEX
#\\begin{adjustbox}{width=10 cm,height= 3cm, keepaspectratio}
	## draw reference line
	my $YSCALE=20;
	if($image{'strand'}=~/\+/) { $image{'genomicStart'}=$image{'start'}; }
	else { $image{'genomicStart'}=$image{'end'}; }
	print OUTFILE "\\draw (0 mm, -2 mm) -- coordinate (x axis mid) (100 mm, -2 mm);\n";
	my $x=0;
	print OUTFILE "\\foreach \\x in {$x,".($x+20).",...,100}\n";
	print OUTFILE "\\draw (\\x mm, -2 mm) -- (\\x mm, -4 mm)\n";
	print OUTFILE "node[anchor=north] {};\n";
	for(my $i=$x; $i<=100; $i+=20) {
		if($i==$x) {
			print OUTFILE "\\draw [black] (".($i-4)." mm, -5 mm) (".($i-4)." mm, -6 mm) node {$image{'chr'}:$image{'genomicStart'}};\n";
		}
		else {
			print OUTFILE "\\draw [black] ($i mm, -5 mm) ($i mm, -6 mm) node {$image{'genomicStart'}};\n";
		}
		if($image{'strand'}=~/\+/) { $image{'genomicStart'} += sprintf("%0.0f", $image{'width'}*(20/100)); }
		else { $image{'genomicStart'} -= sprintf("%0.0f", $image{'width'}*(20/100)); }
	}

	## print information about tick scale
	my $TICKS=sprintf("%0.0f", $image{'width'}*(20/100));
	print OUTFILE "\\draw [black] (50 mm, -10 mm) (50 mm, -11 mm) node {Ticks on x-axis are at $TICKS nucleotides interval};\n";
	

	my $max_y=0;
	## draw read profile (bigWig)
	foreach my $l(@{$track{'bigWig'}}) {
		if($image{'strand'}=~/\+/) {
			my @t=split(/\s+/,$l);
			my $x1=(($t[1]-$image{'start'})/$image{'width'})*100;
			my $x2=(($t[2]-$image{'start'})/$image{'width'})*100;
			my $y=($t[3]/$image{'maxExpr'})*$YSCALE;
			print OUTFILE "\\draw [thick,black,fill=black] ($x1 mm, 0 mm) rectangle ($x2 mm, $y mm);";
			$max_y=$y+5 if($y>$max_y);
		}
		elsif($image{'strand'}=~/\-/) {
			my @t=split(/\s+/,$l);
			my $x1=(($image{'end'}-$t[2])/$image{'width'})*100;
			my $x2=(($image{'end'}-$t[1])/$image{'width'})*100;
			my $y=($t[3]/$image{'maxExpr'})*$YSCALE;
			print OUTFILE "\\draw [thick,black,fill=black] ($x1 mm, 0 mm) rectangle ($x2 mm, $y mm);";
			$max_y=$y+5 if($y>$max_y);
		}
	}

	## draw read profile (wig)
	foreach my $l(@{$track{'wig'}}) {
		if($image{'strand'}=~/\+/) {
			my @t=split(/\s+/,$l);
			my $x1=(($t[0]-$image{'start'})/$image{'width'})*100;
			my $x2=((($t[0]-$image{'start'})+1)/$image{'width'})*100;
			my $y=($t[1]/$image{'maxExpr'})*$YSCALE;
			print OUTFILE "\\draw [thick,black,fill=black] ($x1 mm, 0 mm) rectangle ($x2 mm, $y mm);";
			$max_y=$y+5 if($y>$max_y);
		}
		elsif($image{'strand'}=~/\-/) {
			my @t=split(/\s+/,$l);
			my $x1=(($image{'end'}-$t[0])/$image{'width'})*100;
			my $x2=((($image{'end'}-$t[0])+1)/$image{'width'})*100;
			my $y=($t[1]/$image{'maxExpr'})*$YSCALE;
			print OUTFILE "\\draw [thick,black,fill=black] ($x1 mm, 0 mm) rectangle ($x2 mm, $y mm);";
			$max_y=$y+5 if($y>$max_y);
		}
	}

	## sort blockbuster track
	my %sortBlockbusterTrack=();
	for(my $i=0; $i<scalar(@{$track{'blockbuster'}}); $i++) {
		my @t=split(/\s+/, $track{'blockbuster'}[$i]);
		$sortBlockbusterTrack{$i}=$t[1];
	}	

	my $max_end=0; my $default_y=$max_y; my $block_height=0.5;
	## draw blocks within the read profile
	if($image{'strand'}=~/\+/) {
		#foreach my $l(@{$track{'blockbuster'}}) {
		foreach my $index(sort { $sortBlockbusterTrack{$a} <=> $sortBlockbusterTrack{$b} } keys(%sortBlockbusterTrack)) {
			my $l=$track{'blockbuster'}[$index];
			my @t=split(/\s+/,$l);
			my $x1=(($t[1]-$image{'start'})/$image{'width'})*100;
			my $x2=(($t[2]-$image{'start'})/$image{'width'})*100;
			my $expr=$t[3];
			$expr=~s/^.+\_//g;
			if(defined($sizeFactor)) {
				$expr=sprintf("%0.0f", $expr/$sizeFactor);
			}

			my $y1=(); my $y2=();
			if($x1>$max_end) { $y1=$default_y; $y2=$y1+$block_height; }
			else { $y1=$max_y; $y2=$y1+$block_height; }
			print OUTFILE "\\draw [thick,red,fill=red] ($x1 mm, $y1 mm) rectangle ($x2 mm, $y2 mm);";
			print OUTFILE "\\draw [black] ($x1 mm, $y1 mm) ($x1 mm, ".($y1+2)." mm) node {$expr};\n";
			$max_y=$y2;
			$max_end=$x2;
		}
	}
	elsif($image{'strand'}=~/\-/) {
		#foreach my $l(@{$track{'blockbuster'}}) {
		foreach my $index(reverse sort { $sortBlockbusterTrack{$a} <=> $sortBlockbusterTrack{$b} } keys(%sortBlockbusterTrack)) {
			my $l=$track{'blockbuster'}[$index];
			my @t=split(/\s+/,$l);
			my $x1=(($image{'end'}-$t[2])/$image{'width'})*100;
			my $x2=(($image{'end'}-$t[1])/$image{'width'})*100;
			my $expr=$t[3];
			$expr=~s/^.+\_//g;
			if(defined($sizeFactor)) {
				$expr=sprintf("%0.0f", $expr/$sizeFactor);
			}

			my $y1=(); my $y2=();
			if($x1>$max_end) { $y1=$default_y; $y2=$y1+$block_height; }
			else { $y1=$max_y; $y2=$y1+$block_height; }
			print OUTFILE "\\draw [thick,red,fill=red] ($x1 mm, $y1 mm) rectangle ($x2 mm, $y2 mm);";
			print OUTFILE "\\draw [black] ($x1 mm, $y1 mm) ($x1 mm, ".($y1+2)." mm) node {$expr};\n";
			$max_y=$y2;
			$max_end=$x2;
		}
	}

print OUTFILE <<LATEX;
\\end{tikzpicture}
\\end{adjustbox}
\\end{document}
LATEX
close OUTFILE;

## compile the tex file
system("pdflatex $tmpFilePrefix"."$pdffile.tex > $tmpFilePrefix"."tmp");

## delete intermediate files.
system("rm $tmpFilePrefix"."$pdffile.aux");
system("rm $tmpFilePrefix"."$pdffile.tex");
system("rm $tmpFilePrefix"."$pdffile.log");
system("rm $tmpFilePrefix"."tmp");
print "$tmpFilePrefix$pdffile.pdf\n";

if($png) {
	system("mogrify -trim -format png $tmpFilePrefix$pdffile.pdf");
	system("rm $tmpFilePrefix$pdffile.pdf");
}
$count++;
}

exit;
