#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#########################################################################################
## parse input options
use vars qw($txtFile $htmlFile $maxAlnPerFile $webWorkPath $ucscTrackFile1 $ucscTrackFile2 $webURL $fTitle $sTitle $mode $help);
$webURL="http://rth.dk/resources/dba";
$webWorkPath=".";
$fTitle="First";
$sTitle="Second";
$mode="1";
$maxAlnPerFile=10;

GetOptions ("t=s"  => \$txtFile,
            "l=s"  => \$htmlFile,
						"u=s"  => \$webWorkPath,
						"i=s"  => \$ucscTrackFile1,
						"j=s"  => \$ucscTrackFile2,
            "r=s"  => \$webURL,
            "v=s"  => \$fTitle,
            "w=s"  => \$sTitle,
            "x=i"  => \$maxAlnPerFile,
						"m=i"  => \$mode,
						"help" => \$help,
            "h"    => \$help);

usage() if($help || !$htmlFile);

#########################################################################################
sub usage {
	print STDERR "\nProgram: txt2html.pl (convert alignment result in text format from plotdeepBlockAlign.pl to html format)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: txt2html.pl -t <file> -l <file> [OPTIONS]\n";
	print STDERR " -t <file>     [plotdeepBlockAlign.pl output in text format | STDIN]\n";
	print STDERR " -l <file>     [output HTML file]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -u <string>   [web path to UCSC track and image files (default: null)]\n";
	print STDERR " -i <file>     [UCSC track file corresponding to query block group(s)]\n";
	print STDERR " -j <file>     [UCSC track file corresponding to subject block group(s)]\n";
	print STDERR " -r <string>   [web path to the home (default: http://rth.dk/resources/dba/)]\n";
  print STDERR " -v <string>   [first block group title (default: first)]\n";
	print STDERR " -w <string>   [second block groups title (default: second)]\n";
	print STDERR " -x <int>      [maximum alignments per output html file (default: 10)]\n";
	print STDERR " -m <int>      [1: comparative search; 2: database search]\n";
	print STDERR "               [3: preprocessing comparative search; 4: preprocessing database search]\n";
	print STDERR " -h            [print this helpful message]\n\n";
	exit(-1);
}

#########################################################################################

## customize result path of DBA web server based on jobid
my @t=split(/\//,$webWorkPath);
my $DBA_RESULT_URL="http://rth.dk/resources/dba/results/".$t[scalar(@t)-1];

## open text file
my $INFILE=();
if(defined($txtFile)) {
	if($txtFile=~/\.gz$/) { open($INFILE, "gunzip -c $txtFile |") || die $!; }
	else { open($INFILE, $txtFile) || die $!; }
}
else { $INFILE = *STDIN; }

## read alignment information from text file
my %alnInfo=(); my %nonAlnInfo=(); my $ALIGNMENT_COUNTER=0; my $NON_ALIGNMENT_COUNTER=0; my $TITLE="";
## preprocessing information from text file
my %preprocessingInfo=();

foreach my $l(<$INFILE>) {
	chomp($l);
	my @t=split(/\t+/, $l);
	if($l=~/^\#title\:/) {
		$TITLE.="$t[1] ";
	}
	elsif($l=~/^DBA_f1/) { $preprocessingInfo{'DBA_f1'}{$t[1]}=$t[2]; }
	elsif($l=~/^DBA_f2/) { $preprocessingInfo{'DBA_f2'}{$t[1]}=$t[2]; }
	elsif($l=~/^\#alignment information/) {
		$ALIGNMENT_COUNTER++;
	}
	elsif($l=~/first block group\:/) {
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'name'} = $t[1];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'chr'} = $t[2];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'start'} = $t[3];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'end'} = $t[4];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'strand'} = $t[5];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'reads'} = $t[6];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'tags'} = $t[7];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'blocks'} = $t[8];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'entropy'} = $t[9];
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'anno'} = $t[10];
		if($t[11]=~/^\.$/) { $alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'loci'} = "intergene"; }
		else { $alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'loci'} = $t[11]; }
		$alnInfo{$ALIGNMENT_COUNTER}{'fBG'}{'uniqReads'} = $t[12];
	}
	elsif($l=~/second block group\:/) {
		$t[1]=~s/\_cluster.+\(/ (/ig;
		$t[1]=~s/\_bg.+\(/ (/ig;
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'name'} = $t[1];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'chr'} = $t[2];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'start'} = $t[3];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'end'} = $t[4];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'strand'} = $t[5];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'reads'} = $t[6];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'tags'} = $t[7];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'blocks'} = $t[8];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'entropy'} = $t[9];
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'anno'} = $t[10];
		if($t[11]=~/^\.$/) { $alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'loci'} = "intergene"; }
		else { $alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'loci'} = $t[11]; }
		$alnInfo{$ALIGNMENT_COUNTER}{'sBG'}{'uniqReads'} = $t[12];
	}
	elsif($l=~/alignment\:/) {
		$alnInfo{$ALIGNMENT_COUNTER}{'score'} = $t[1];
		$alnInfo{$ALIGNMENT_COUNTER}{'bgImage'} = $t[2];
		if($mode==1) {
			$alnInfo{$ALIGNMENT_COUNTER}{'pvalue'} = $t[3];
		}
		elsif($mode==2) {
			foreach(@t[3..(scalar(@t)-1)]) { push(@{$alnInfo{$ALIGNMENT_COUNTER}{'bImage'}}, $_);	}
		}
	}
	elsif($l=~/^\#non\-aligned\:/) {
		$NON_ALIGNMENT_COUNTER++;
		$nonAlnInfo{$NON_ALIGNMENT_COUNTER}{'fBG'}{'name'}=$t[1];
		$nonAlnInfo{$NON_ALIGNMENT_COUNTER}{'fBG'}{'coor'}=$t[2];
  }
}

## format alignment information as HTML output
if($mode==1) {

my $HTMLFILE=(); my %seen=();
my $COUNTER=(); my $HTML_FILE_COUNTER=0;
my $START_COUNTER=(); my $END_COUNTER=();

## create result file(s) in HTML format
for($COUNTER=1; $COUNTER<=$ALIGNMENT_COUNTER; $COUNTER++) {

## print header, if atleast one significant alignment is observed
if(scalar(keys(%alnInfo))>=1 && ($COUNTER%$maxAlnPerFile)==1) {

$HTML_FILE_COUNTER++; $START_COUNTER=$COUNTER;
close $HTMLFILE if(defined($HTMLFILE));
open($HTMLFILE, ">$htmlFile$HTML_FILE_COUNTER.tpl") || die $!; 

print $HTMLFILE "<p>$TITLE</p>\n";

## print HTML link to browse over result files
for(my $i=1; $i<=sprintf("%0.0f",$ALIGNMENT_COUNTER/$maxAlnPerFile); $i++) {
	if($i==$HTML_FILE_COUNTER) {
	  print $HTMLFILE "<a href=\"$DBA_RESULT_URL/$i\"><b>$i</b></a>&nbsp;";
	}
	else {
	  print $HTMLFILE "<a href=\"$DBA_RESULT_URL/$i\">$i</a>&nbsp;";
	}
}

print $HTMLFILE <<HTML
<table id=resultTable class=resultTable border=0>
	<tr class="header"><th>S.No.</th><th>Genomic loci</th><th>$fTitle Sample<a target="_blank" href="$webURL/help#fBG">[?]</a></th><th>$sTitle Sample<a target="_blank" href="$webURL/help#sBG">[?]</a></th><th>Score<a target="_blank" href="$webURL/help#score">[?]</a></th><th>Description<a target="_blank" href="$webURL/help#description">[?]</a></th><th>Read profile<a target="_blank" href="$webURL/help#readProfile">[?]</a></th><tr>
HTML

}

print $HTMLFILE <<HTML;
	<tr class="result"><td>$COUNTER</td><td><a name="coordinate" style="cursor:pointer;" onClick='managePlot("readProfile$COUNTER");'><font color="#305080"><u>$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}</u></font></a></td><td>$alnInfo{$COUNTER}{'fBG'}{'name'}</td><td>$alnInfo{$COUNTER}{'sBG'}{'name'}</td><td>$alnInfo{$COUNTER}{'score'}</td>
		<td><img onClick='managePlot("description$COUNTER");' id="img_description$COUNTER" src="http://rth.dk/resources/dba/img/add.png"></td>
		<td><img onClick='managePlot("readProfile$COUNTER");' id="img_readProfile$COUNTER" src="http://rth.dk/resources/dba/img/add.png"></td>
	</tr>
	<tr><td colspan=\"7\">
		<table id="description$COUNTER" class="detail">
			<tr class="sub_header"><th class="emptyth" rowspan="10"></th><th colspan="2">$fTitle block group</th><th colspan="2">$sTitle block group</th></tr>
			<tr><th>name<a target="_blank" href="$webURL/help#name">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'name'}</td><th>name<a target="_blank" href="$webURL/help#name">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'name'}</td></tr>
HTML
			if(defined($webWorkPath) && defined($ucscTrackFile1) && defined($ucscTrackFile2)) {
print $HTMLFILE <<HTML;
			<tr><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org=hg19&position=$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}&hgt.customText=$webWorkPath/$ucscTrackFile1" target="_blank">$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}</a>($alnInfo{$COUNTER}{'fBG'}{'strand'})</td><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org=hg19&position=$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}&hgt.customText=$webWorkPath/$ucscTrackFile2" target="_blank">$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}</a>($alnInfo{$COUNTER}{'sBG'}{'strand'})</td></tr>
HTML
			}
			else {
print $HTMLFILE <<HTML;
			<tr><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}($alnInfo{$COUNTER}{'fBG'}{'strand'})</td><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}($alnInfo{$COUNTER}{'sBG'}{'strand'})</td></tr>
HTML
			}
print $HTMLFILE <<HTML;
			<tr><th>reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'reads'}</td><th>reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'reads'}</td></tr>
			<tr><th>tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'tags'}</td><th>tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'tags'}</td></tr>
			<tr><th>blocks<a target="_blank" href="$webURL/help#blocks">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'blocks'}</td><th>blocks<a target="_blank" href="$webURL/help#blocks">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'blocks'}</td></tr>
			<tr><th>entropy<a target="_blank" href="$webURL/help#entropy">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'entropy'}</td><th>entropy<a target="_blank" href="$webURL/help#entropy">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'entropy'}</td></tr>
			<tr><th>ncRNA<a target="_blank" href="$webURL/help#ncrna">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'anno'}</td><th>ncRNA<a target="_blank" href="$webURL/help#ncrna">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'anno'}</td></tr>
			<tr><th>loci<a target="_blank" href="$webURL/help#loci">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'loci'}</td><th>loci<a target="_blank" href="$webURL/help#loci">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'loci'}</td></tr>
			<tr><th>unique reads (%)<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'uniqReads'}</td><th>unique reads (%)<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'uniqReads'}</td></tr>
		</table>
		<table id="readProfile$COUNTER" class="detail">
			<tr class="sub_header"><th rowspan="2"></th><th>Read profile<a target="_blank" href="$webURL/help#readProfile">[?]</a></th></tr>
			<tr><td>
				<div id="align">
HTML
				print $HTMLFILE "<img id=\"img_readProfile\" src=\"$webWorkPath/$alnInfo{$COUNTER}{'bgImage'}\" border=0>";
print $HTMLFILE <<HTML;
				</div>
			</td></tr>
		</table>
	</td></tr>
HTML

if(($COUNTER%$maxAlnPerFile)==0 || $COUNTER==$ALIGNMENT_COUNTER) {
$END_COUNTER=$COUNTER;
print $HTMLFILE <<HTML;
	<tr><td colspan="7"><input type=hidden id="ALIGNMENT_COUNTER_START" value=$START_COUNTER></td></tr>
	<tr><td colspan="7"><input type=hidden id="ALIGNMENT_COUNTER_END" value=$END_COUNTER></td></tr>
</table>
HTML
}
}

close $HTMLFILE if(defined($HTMLFILE));
if($ALIGNMENT_COUNTER==0) {
	open($HTMLFILE, ">$htmlFile"."1.tpl") || die $!; 
	print $HTMLFILE "<p>$TITLE</p>\n";
	close $HTMLFILE;
}
}
elsif($mode==2) {

my $HTMLFILE=(); my %seen=();
my $COUNTER=(); my $BG_COUNTER=(); my $HTML_FILE_COUNTER=0;
my $START_COUNTER=(); my $END_COUNTER=();

## create result file(s) in HTML format
for($COUNTER=1; $COUNTER<=$ALIGNMENT_COUNTER; $COUNTER++) {

## print header, if atleast one significant alignment is observed
if(scalar(keys(%alnInfo))>=1 && ($COUNTER%$maxAlnPerFile)==1) {

$HTML_FILE_COUNTER++; $START_COUNTER=$COUNTER;
close $HTMLFILE if(defined($HTMLFILE));
open($HTMLFILE, ">$htmlFile$HTML_FILE_COUNTER.tpl") || die $!;

print $HTMLFILE "<p>$TITLE</p>\n";

## print HTML link to browse over result files
for(my $i=1; $i<=sprintf("%0.0f",$ALIGNMENT_COUNTER/$maxAlnPerFile); $i++) {
	if($i==$HTML_FILE_COUNTER) { 
		print $HTMLFILE "<a href=\"$DBA_RESULT_URL/$i\"><b>$i</b></a>&nbsp;";
	}
	else {
		print $HTMLFILE "<a href=\"$DBA_RESULT_URL/$i\">$i</a>&nbsp;";
	}
}

print $HTMLFILE <<HTML;
<table id=resultTable class=resultTable border=0>
	<tr class="header"><th>S.No.</th><th>$fTitle block group<a target="_blank" href="$webURL/help#fBG">[?]</a></th><th>$sTitle block group<a target="_blank" href="$webURL/help#sBG">[?]</a></th><th>Score<a target="_blank" href="$webURL/help#score">[?]</a></th><th>Description<a target="_blank" href="$webURL/help#description">[?]</a></th><th>Block group alignment<a target="_blank" href="$webURL/help#bgAlign">[?]</a></th><th>Block alignment<a target="_blank" href="$webURL/help#bAlign">[?]</a></th><tr>
HTML

}

if(!defined($seen{"$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}"})) {
	$BG_COUNTER++;
	print $HTMLFILE "<tr class=\"result\"><td>$BG_COUNTER</td><td>$alnInfo{$COUNTER}{'fBG'}{'name'}</td><td>$alnInfo{$COUNTER}{'sBG'}{'name'}</td><td>$alnInfo{$COUNTER}{'score'}</td>\n";
	$seen{"$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}"}=1;
}
else {
	print $HTMLFILE "<tr class=\"result\"><td>&nbsp;</td><td></td><td>$alnInfo{$COUNTER}{'sBG'}{'name'}</td><td>$alnInfo{$COUNTER}{'score'}</td>\n";
}

print $HTMLFILE <<HTML;
		<td><img onClick='managePlot("description$COUNTER");' id="img_description$COUNTER" src="http://rth.dk/resources/dba/img/add.png"></td>
		<td><img onClick='managePlot("bg_align$COUNTER");' id="img_bg_align$COUNTER" src="http://rth.dk/resources/dba/img/add.png"></td>
		<td><img onClick='managePlot("b_align$COUNTER");' id="img_b_align$COUNTER" src="http://rth.dk/resources/dba/img/add.png"></td>
	</tr>
	<tr><td colspan=\"7\">
		<table id="description$COUNTER" class="detail">
			<tr class="sub_header"><th class="emptyth" rowspan="10"></th><th colspan="2">$fTitle block group</th><th colspan="2">$sTitle block group</th></tr>
			<tr><th>name<a target="_blank" href="$webURL/help#name">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'name'}</td><th>name<a target="_blank" href="$webURL/help#name">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'name'}</td></tr>
HTML

			if(defined($webWorkPath) && defined($ucscTrackFile1) && defined($ucscTrackFile2)) {
print $HTMLFILE <<HTML;
			<tr><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org=hg19&position=$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}&hgt.customText=$webWorkPath/$ucscTrackFile1" target="_blank">$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}</a>($alnInfo{$COUNTER}{'fBG'}{'strand'})</td><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org=hg19&position=$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}&hgt.customText=$webWorkPath/$ucscTrackFile2" target="_blank">$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}</a>($alnInfo{$COUNTER}{'sBG'}{'strand'})</td></tr>
HTML
			}
			else {
print $HTMLFILE <<HTML;
			<tr><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'chr'}:$alnInfo{$COUNTER}{'fBG'}{'start'}-$alnInfo{$COUNTER}{'fBG'}{'end'}($alnInfo{$COUNTER}{'fBG'}{'strand'})</td><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'chr'}:$alnInfo{$COUNTER}{'sBG'}{'start'}-$alnInfo{$COUNTER}{'sBG'}{'end'}($alnInfo{$COUNTER}{'sBG'}{'strand'})</td></tr>
HTML
			}

print $HTMLFILE <<HTML;
			<tr><th>reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'reads'}</td><th>reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'reads'}</td></tr>
			<tr><th>tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'tags'}</td><th>tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'tags'}</td></tr>
			<tr><th>blocks<a target="_blank" href="$webURL/help#blocks">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'blocks'}</td><th>blocks<a target="_blank" href="$webURL/help#blocks">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'blocks'}</td></tr>
			<tr><th>entropy<a target="_blank" href="$webURL/help#entropy">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'entropy'}</td><th>entropy<a target="_blank" href="$webURL/help#entropy">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'entropy'}</td></tr>
			<tr><th>ncRNA<a target="_blank" href="$webURL/help#ncrna">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'anno'}</td><th>ncRNA<a target="_blank" href="$webURL/help#ncrna">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'anno'}</td></tr>
			<tr><th>loci<a target="_blank" href="$webURL/help#loci">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'loci'}</td><th>loci<a target="_blank" href="$webURL/help#loci">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'loci'}</td></tr>
			<tr><th>unique reads (%)<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$alnInfo{$COUNTER}{'fBG'}{'uniqReads'}</td><th>unique reads (%)<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$alnInfo{$COUNTER}{'sBG'}{'uniqReads'}</td></tr>
		</table>
		<table id="bg_align$COUNTER" class="detail">
			<tr class="sub_header"><th rowspan="2"></th><th>Block group alignment<a target="_blank" href="$webURL/help#bgAlign">[?]</a></th></tr>
			<tr><td>
				<div id="align">
					<img id="img_bg_align" src="$webWorkPath/$alnInfo{$COUNTER}{'bgImage'}" border=0>
				</div>
			</td></tr>
		</table>
		<table id="b_align$COUNTER" class="detail">
			<tr class="sub_header"><th rowspan="10"></th><th>Block alignment<a target="_blank" href="$webURL/help#bAlign">[?]</a></th></tr>
HTML

			foreach(my $i=0; $i<scalar(@{$alnInfo{$COUNTER}{'bImage'}}); $i++) {
				if($i%2==0) {
print $HTMLFILE <<HTML;
			<tr><td>
				<div id="align">
HTML
				}
print $HTMLFILE <<HTML;
					<img id="img_b_align" src="$webWorkPath/$alnInfo{$COUNTER}{'bImage'}[$i]" border=0>
HTML
				if($i%2!=0 || $i==(scalar(@{$alnInfo{$COUNTER}{'bImage'}})-1)) {
print $HTMLFILE <<HTML;
				</div>
			</td></tr>
HTML
				}
			}
print $HTMLFILE <<HTML;
		</table>
	</td></tr>
HTML


if(($COUNTER%$maxAlnPerFile)==0 || $COUNTER==$ALIGNMENT_COUNTER) {
$END_COUNTER=$COUNTER;
print $HTMLFILE <<HTML;
	<tr><td colspan="7"><input type=hidden id="ALIGNMENT_COUNTER_START" value=$START_COUNTER></td></tr>
	<tr><td colspan="7"><input type=hidden id="ALIGNMENT_COUNTER_END" value=$END_COUNTER></td></tr>
</table>
HTML
}

if(keys(%nonAlnInfo) && ($COUNTER%$maxAlnPerFile)==0) {
print $HTMLFILE "<p>No significant hit is observed for ".(keys(%nonAlnInfo))." block group(s)</p>\n";
print $HTMLFILE <<HTML;
<table id=resultTable class=resultTable border=0>
	<tr class="header"><th>$fTitle block group<a target="_blank" href="$webURL/help#fBG">[?]</a></th><th>coordinate<a target="_blank" href="$webURL/help#coordinate">[?]</a></th><tr>
HTML

	foreach $NON_ALIGNMENT_COUNTER(sort { $a <=> $b } keys(%nonAlnInfo)) {
print $HTMLFILE <<HTML;
	<tr class="result"><td>$nonAlnInfo{$NON_ALIGNMENT_COUNTER}{'fBG'}{'name'}</td><td>$nonAlnInfo{$NON_ALIGNMENT_COUNTER}{'fBG'}{'coor'}</td></tr>
HTML
	}
print $HTMLFILE <<HTML;
</table>
HTML
}

}

close $HTMLFILE if(defined($HTMLFILE));
if($ALIGNMENT_COUNTER==0) {
	open($HTMLFILE, ">$htmlFile"."1.tpl") || die $!; 
	print $HTMLFILE "<p>$TITLE</p>\n";
	close $HTMLFILE;
}
}
elsif($mode==3) {

my $HTMLFILE=();
open($HTMLFILE, ">$htmlFile.tpl") || die $!;

print $HTMLFILE <<HTML;
<table id=outerTable class=outerTable border=0>
	<tr><td valign="top">
		<table class=preprocessingTable border=0>
			<tr class="header"><th colspan="2">Input BED file</th></tr>
			<tr><th>Reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'reads'}</td></tr>
			<tr><th>Uniquely mapped reads<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'uniqueMappingReads'}</td></tr>
			<tr><th>Tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'tags'}</td></tr>
			<tr><th>Unqiuely mapped tags<a target="_blank" href="$webURL/help#uniqTags">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'uniqueMappingTags'}</td></tr>
		</table>
	</td>
	<td>
		<table class=preprocessingTable border=0>
			<tr class="header"><th colspan="2">Block groups processed from input BED file<a target="_blank" href="$webURL/help#preprocessing">[?]</a></th></tr>
			<tr><th>miRNA<a target="_blank" href="$webURL/help#miRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'miRNA'}</td></tr>
			<tr><th>snoRNA<a target="_blank" href="$webURL/help#snoRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'snoRNA'}</td></tr>
			<tr><th>tRNA<a target="_blank" href="$webURL/help#tRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'tRNA'}</td></tr>
			<tr><th>Other ncRNA<a target="_blank" href="$webURL/help#other">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'other ncRNA'}</td></tr>
			<tr><th>Unannotated<a target="_blank" href="$webURL/help#unannotated">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'Unannotated'}</td></tr>
			<tr><th>Total<a target="_blank" href="$webURL/help#total">[?]</a></th><td><a href="$webWorkPath/DBA_f1.map.clusters.flagged.gz">$preprocessingInfo{'DBA_f1'}{'Total'}</a></td></tr>
		</table>
	</td></tr>
</table>
HTML
close $HTMLFILE;
}
elsif($mode==4) {

my $HTMLFILE=();
open($HTMLFILE, ">$htmlFile.tpl") || die $!;

print $HTMLFILE <<HTML;
<table id=outerTable class=outerTable border=0>
	<tr><td valign="top">
		<table class=preprocessingTable border=0>
			<tr class="header"><th colspan="3">Input BED file</th></tr>
			<tr><th>&nbsp;</th><td>Sample 1</td><td>Sample 2</td></tr>
			<tr><th>Reads<a target="_blank" href="$webURL/help#reads">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'reads'}</td><td>$preprocessingInfo{'DBA_f2'}{'reads'}</td></tr>
			<tr><th>Uniquely mapped reads<a target="_blank" href="$webURL/help#uniqReads">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'uniqueMappingReads'}</td><td>$preprocessingInfo{'DBA_f2'}{'uniqueMappingReads'}</td></tr>
			<tr><th>Tags<a target="_blank" href="$webURL/help#tags">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'tags'}</td><td>$preprocessingInfo{'DBA_f2'}{'tags'}</td></tr>
			<tr><th>Unqiuely mapped tags<a target="_blank" href="$webURL/help#uniqTags">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'uniqueMappingTags'}</td><td>$preprocessingInfo{'DBA_f2'}{'uniqueMappingTags'}</td></tr>
		</table>
	</td>
	<td>
		<table class=preprocessingTable border=0>
			<tr class="header"><th colspan="3">Block groups processed from input BED file<a target="_blank" href="$webURL/help#preprocessing">[?]</a></th></tr>
			<tr><th>&nbsp;</th><td>Sample 1</td><td>Sample 2</td></tr>
			<tr><th>miRNA<a target="_blank" href="$webURL/help#miRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'miRNA'}</td><td>$preprocessingInfo{'DBA_f2'}{'miRNA'}</td></tr>
			<tr><th>snoRNA<a target="_blank" href="$webURL/help#snoRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'snoRNA'}</td><td>$preprocessingInfo{'DBA_f2'}{'snoRNA'}</td></tr>
			<tr><th>tRNA<a target="_blank" href="$webURL/help#tRNA">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'tRNA'}</td><td>$preprocessingInfo{'DBA_f2'}{'tRNA'}</td></tr>
			<tr><th>Other ncRNA<a target="_blank" href="$webURL/help#other">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'other ncRNA'}</td><td>$preprocessingInfo{'DBA_f2'}{'other ncRNA'}</td></tr>
			<tr><th>Unannotated<a target="_blank" href="$webURL/help#unannotated">[?]</a></th><td>$preprocessingInfo{'DBA_f1'}{'Unannotated'}</td><td>$preprocessingInfo{'DBA_f2'}{'Unannotated'}</td></tr>
			<tr><th>Total<a target="_blank" href="$webURL/help#total">[?]</a></th><td><a href="$webWorkPath/DBA_f1.map.clusters.flagged.gz">$preprocessingInfo{'DBA_f1'}{'Total'}</a></td><td><a href="$webWorkPath/DBA_f2.map.clusters.flagged.gz">$preprocessingInfo{'DBA_f2'}{'Total'}</a></td></tr>
		</table>
	</td></tr>
</table>
HTML

close $HTMLFILE;
}

exit;
