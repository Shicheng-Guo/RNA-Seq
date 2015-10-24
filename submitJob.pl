#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

##########################################################################################################
# parse input options
use vars qw($qUniqkeyword $qCommand $qSystem $help);
$qSystem="dcsc";

GetOptions ("k=s"  => \$qUniqkeyword,
            "c=s"  => \$qCommand,
            "q=s"  => \$qSystem,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$qUniqkeyword || !$qCommand);

##########################################################################################################
sub usage {
	print STDERR "\nProgram: submitJob.pl (submit job to queue)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: submitJob.pl -k <string> -c <string>\n";
	print STDERR " -k <string>  [unique keyword to associate with each script submitted in queue]\n";
	print STDERR " -c <string>  [command to execute within double quotes]\n";
	print STDERR " -q <string>  [queue system (rth or dcsc; default: dcsc)]\n";
	print STDERR " -h           [help]\n\n";
	exit(-1);
}
##########################################################################################################

open(OUTFILE, ">theScript_$qUniqkeyword.sh");

if($qSystem=~/dcsc/) {
	print OUTFILE "# @ output = theScript_$qUniqkeyword.out\n";
	print OUTFILE "# @ error = theScript_$qUniqkeyword.err\n";
	print OUTFILE "# @ job_type = serial\n";
	print OUTFILE "# @ resources   = ConsumableCpus(1) ConsumableMemory(12 gb)\n";
	print OUTFILE "# @ env_copy = all\n";
	print OUTFILE "# @ class = rth\n";
	print OUTFILE "# @ user_priority = 100\n";
	print OUTFILE "# @ queue\n";
	print OUTFILE "#memory limit in kb\n";
	print OUTFILE "ulimit -v 12000000\n";
	print OUTFILE "#scratch directory\n";
	print OUTFILE "cd /scratch/\$LOADL_STEP_ID\n";
	print OUTFILE "#directory from where the job was submitted\n";
	print OUTFILE "cd \$LOADL_STEP_INITDIR\n";
	$qCommand=~s/\"//g;
	print OUTFILE "$qCommand\n";
	close OUTFILE;

	## submit job to queue
	system("llsubmit theScript_$qUniqkeyword.sh");
}
else {
	print OUTFILE "#!/bin/bash\n";
	print OUTFILE "#PBS -l nodes=1:ppn=4\n";
	$qCommand=~s/\"//g;
	print OUTFILE "$qCommand\n";
	close OUTFILE;

	## submit job to queue
	system("qsub theScript_$qUniqkeyword.sh");
}

exit;
