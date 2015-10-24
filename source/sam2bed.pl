#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Getopt::Std;
use List::Util;
use Cwd;
use IO::File;
use POSIX qw(tmpnam);

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $inputFile $outputFile $outFormatArf $bam);


# -----------------------------------------------------------------------------
# OPTIONS

GetOptions (
"i=s"       => \$inputFile,
"o=s"       => \$outputFile,
"f=s"       => \$outFormatArf,
"b"         => \$bam,
"help"      => \$help,
"h"         => \$help);
usage() if ($help || !$inputFile || !$outputFile);


# -----------------------------------------------------------------------------
# MAIN

printHeader();

sam2bed($inputFile, $outputFile);



# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
    print STDERR "\nusage: sam2bed.pl -i <file> -o <file> [-b]\n";
    print STDERR "bring mapping output to bed format\n";
    print STDERR "\n";
    print STDERR "[INPUT]\n";
    print STDERR " -i <file>    mapped reads file\n";
    print STDERR " -b           bam file\n";
		print STDERR " -f           output format is arf (mirdeep2)\n";
    print STDERR " -o <file>    output file\n";
    print STDERR " -h <file>    this (usefull) help message\n";
    print STDERR "[VERSION]\n";
    print STDERR " 02-22-2012\n";
    print STDERR "[BUGS]\n";
    print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
    print STDERR "\n";
    exit(-1);
}

sub prettyTime{
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month,
    $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    return "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}

sub printHeader{
    print "# sam2bed.pl started " . prettyTime() . "\n";
    print "# input file: $inputFile\n";
    print "#\n";
}

sub sam2bed{
    my ($inputFile, $outputFile) = @_;
    my %tags = (); my $c = 0; my %coor = ();
    if(!$bam){
        open(FILE, "<$inputFile") || die "cannot open $inputFile\n";
    }elsif($inputFile=~/gz$/) {
        open(FILE, "zless $inputFile | samtools view - | ") || die "cannot open $inputFile\n";
    }else{
        open(FILE, "samtools view $inputFile | ") || die "cannot open $inputFile\n";
    }

    while(<FILE>){
        chomp;
        my @entry = split(/\t/,$_);
        next if ($_ =~ /^[@]/ || $_ =~ /^\s*$/ || $entry[2] eq "*");
        checkFormat($_, "SAM");
        my ($id, $chr, $strand, $seq, $start, $end, $cigar) = sam($_);

        # name tag
        if(!exists($tags{$seq}{tag})){
            $c++;
            my $tagId = $inputFile;
            $tagId=~s/^.+\///g;
            $tagId=~s/\.[bs]am*//g;
            $tags{$seq}{tag} = $tagId."_".$c;
        }
        # count reads
        # the count of distict ids for a seq in bam file is the expression of a read
        #if(!exists($tags{$seq}{reads}{$id})){
            $tags{$seq}{reads}{$id}++;
        #}
        # count loci
        if($cigar !~ /N/){
            $tags{$seq}{loci}{"$chr|$start|$end|$strand"}++;
            $coor{$seq}{$id}{"$chr|$start|$end|$strand"}++;
        }else{
            my @parts = getSpliceParts($cigar, $start);
            foreach my $part (@parts){
                ($start, $end) = split(/\|/, $part);
								if(($end-$start)>40) { print STDERR $_."\n"; }
                $tags{$seq}{loci}{"$chr|$start|$end|$strand"}++;
                $coor{$seq}{$id}{"$chr|$start|$end|$strand"}++;
            }
        }

        foreach my $col (@entry){
            if($col =~ /XA:Z:/ && $col !~ /XA:Z:Q/){
                $col =~ s/XA:Z://;
                my @multi = split(/\;/,$col);
                foreach my $mm (@multi){
                    ($chr, $start, $cigar) = split(/\,/, $mm);
                    if($start =~ /\+/){$strand = "+";}else{$strand = "-";}
                    $start =~ s/[\+\-]//;
                    my $length = 0; my $tmp = $cigar; $tmp =~ s/(\d+)[MD]/$length+=$1/eg;
                    my $end = $start + $length;
                    $start--; $end--;
                    if($cigar !~ /N/){
                        $tags{$seq}{loci}{"$chr|$start|$end|$strand"}++;
                    }
                    else{
                        my @parts = getSpliceParts($cigar, $start);
                        foreach my $part (@parts){
                            ($start, $end) = split(/\|/, $part);
														if(($end-$start)>40) { print STDERR $_."\n"; }
                            $tags{$seq}{loci}{"$chr|$start|$end|$strand"}++;
                        }
                    }
                }
            }
        }
    }
    close(FILE);

		## For genomic loci ($locus) where a sequence ($seq) is mapped,
		##   sequence is uniquely mapped, if it has only one tag id ($id) associated with it
		##   sequence is multi-mapped, if it has more than one tag id ($id) associated with it
    my %chroms = (); my $uid=1;
    foreach my $seq (keys %tags){
        #print "$seq\t".keys(%{$tags{$seq}{loci}})."\t".keys(%{$tags{$seq}{reads}})."\t";
        #my $sum=0;
        #foreach(keys(%{$tags{$seq}{reads}})) {$sum+=$tags{$seq}{reads}{$_};}
        #print "$sum\n";
        foreach my $locus (keys %{$tags{$seq}{loci}}){
            my ($chr, $start, $end, $strand) = split(/\|/,$locus);
            if(!exists($chroms{$chr})){($chroms{$chr}{fileName},$chroms{$chr}{fileHandle}) = openBin();}
            my $file = $chroms{$chr}{fileHandle};
            my $uniqTags=0; my $duplTags=0;
            foreach my $id(keys(%{$tags{$seq}{reads}})) {
                if($tags{$seq}{reads}{$id}==1) {
                    foreach my $coor(keys(%{$coor{$seq}{$id}})) {
                        if($coor=~/^\Q$locus\E$/) { $uniqTags++; }
                    }
                }
                else {
                    foreach my $coor(keys(%{$coor{$seq}{$id}})) {
                        if($coor=~/^\Q$locus\E$/) { $duplTags++; }
                    }
                }
            }
						if($outFormatArf) {
	            if($uniqTags>0) {
									my $ref=`coor2seq.pl -c $chr:$start-$end -s $strand`;
									chomp($ref);
	                print $file "$tags{$seq}{tag}_".$uid."_x$duplTags\t".(length($seq))."\t1\t".(length($seq))."\t$seq\t$chr\t".($end-$start)."\t$start\t$end\t$ref\t$strand\t";
	            }
	            if($duplTags>0) {
	            }
						}
						else {
	            if($uniqTags>0) {
	                print $file "$chr\t$start\t$end\t$tags{$seq}{tag}_1\t".$uniqTags."\t$strand\n";
	            }
	            if($duplTags>0) {
	                print $file "$chr\t$start\t$end\t$tags{$seq}{tag}_".(keys(%{$tags{$seq}{loci}}))."\t".$duplTags / (keys(%{$tags{$seq}{loci}}))."\t$strand\n";
	            }
						}
						$uid++;
        }
    }
    closeBins(%chroms);
    mergeAndSort(%chroms);
    deleteBins(%chroms);
}

sub checkFormat{
    my ($entry, $format) = @_;
    if($format eq "SAM"){
        my @l = split(/\t+/, $entry);
        if($l[3] !~ /\d+/){die "file is not in SAM format -> $entry\n";}
        if($l[1] != 0 && $l[1] != 16 && $l[1] != 4){die "It seems your file consists of paired end reads.\nUnfortunately, we do not support these experiments.\n";}
    }
}

sub sam{
    my ($line) = @_;
    my @entry = split(/\t/,$line);
    my $id = $entry[0];
    my $chr = $entry[2];
    my $strand = $entry[1]?"-":"+";
    my $seq = $entry[9]; if($strand eq "-"){$seq = reverseComplement($seq);}
    my $start = $entry[3];
    my $cigar = $entry[5];
    my $length = 0; my $tmp = $cigar; $tmp =~ s/(\d+)[MD]/$length+=$1/eg;
    my $end = $entry[3] + $length;
    $start--; $end--;
		if(($end-$start)>40) { print STDERR $line."\n"; }

    return ($id, $chr, $strand, $seq, $start, $end, $cigar);
}

sub getSpliceParts{
    my ($cigar, $start) = @_;
    my @parts = ();
    my $length = 0; my $value = "";
    for(my $i = 0; $i <= length($cigar); $i++){
        my $v = substr($cigar,$i,1);
        if($v=~/\d/){$value.=$v;}
        if($v=~/\D/){
            if($v eq "M" || $v eq "D"){
							$length+=$value;$value="";
						}
            elsif($v eq "N"){
                my $end += $start + $length + 1;
                push(@parts, "$start|$end");
                $start = $start + $length + $value;
                $length = 0;$value="";
						}
						else {$value="";}
				}
		}
    my $end += $start + $length + 1;
    push(@parts, "$start|$end");
    return @parts;
}

sub reverseComplement{
    my ($seq) = @_;
    my $revcomp = reverse($seq);
    $revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;
    return $revcomp;
}

sub openBin{
    my ($fileName, $file) = ();
    do {$fileName = tmpnam()} until $file = IO::File->new($fileName, O_RDWR|O_CREAT|O_EXCL);
    return ($fileName, $file);
}

sub deleteBins{
    my (%bins) = @_;
    foreach my $id (keys %bins){
        unlink($bins{$id}{fileName});
    }
}

sub closeBins{
    my (%bins) = @_;
    foreach my $id (keys %bins){
        close($bins{$id}{fileHandle});
    }
}

sub mergeAndSort{
    my (%chroms) = @_;
    my %ids = ();
    foreach my $c (keys %chroms){
        my ($fileName,$fileHandle) = openBin();
        system("nohup sort -k 1,1 -k 6,6 -k 2n,2 -k 3n,3 -o $fileName $chroms{$c}{fileName}  2>&1 1>/dev/null &");
        $ids{$fileName}{fileName} = $fileName;
        $ids{$fileName}{fileHandle} = $fileHandle;
    }
    while(checkFiles(%ids) != 1){
        #print "sleep(10)\t".prettyTime()."\n\n";
        sleep(10);
    }
    closeBins(%ids);

    if(-e $outputFile){system("rm $outputFile");}
    foreach my $id (keys %ids){
        system("cat $ids{$id}{fileName} >> $outputFile");
    }
    deleteBins(%ids);
}

sub checkFiles{
    my (%ids) = @_;
    my $ps = `ps -f`;
    #print "$ps\n";
    foreach my $id (keys %ids){
        if(defined($ps) && $ps =~ /$id/){return 0;};
    }

    return 1;
}

