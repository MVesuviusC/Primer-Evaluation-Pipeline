#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 1/16/20
# Last modified: 
# Title: batchQsub.pl
# Purpose: 
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $input;
my $outDir = "evalOut";
my $proc = 30;
my $blastLoc = "/usr/local/packages/ncbi-blast+-2.7.1/bin/blastn";
my $blastDb = "/local/projects-t3/SerreDLab-3/databases/blast/nt";
my $taxonomyDb = "/local/projects-t3/SerreDLab-3/databases/taxonomy.db";


# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
      	    "input=s"           => \$input,
      	    "outDir=s"          => \$outDir,
	    "proc=i"            => \$proc,
	    "blastLoc=s"        => \$blastLoc,
	    "blastDb=s"         => \$blastDb,
	    "taxDb=s"           => \$taxonomyDb
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);

##############################
# Global variables
##############################

##############################
# Code
##############################
$outDir .= '/';

##############################
### Make header for shell script 

my $qsubHeader = 
    '#!/bin/bash' . "\n" .
    '#$ -cwd' . "\n" . 
    '#$ -P dserre-lab' . "\n" .
    '#$ -l mem_free=40G' . "\n" .
    '#$ -q threaded.q' . "\n" .
    '#$ -pe thread ' . $proc . "\n" .
    '#$ -v PATH';


##############################
### Make base output directory

if(-w "." eq "") { # check for write permissions
      print STDERR "No write permissions to this directory!\n\n";
      die;
}

mkdir $outDir;

if($verbose) {
    print STDERR "Output directory: ", $outDir, "\n\n";
}


##############################
### Make separate shell scripts for each primer set and submit them

open my $inputFH, "$input" or die "Could not open primer input\nWell, crap\n";

while (my $input = <$inputFH>){
    chomp $input;
    my ($primerFName, $primerF, $primerRName, $primerR, $taxa, $level) = split "\t", $input;
    
    # Make output directory for each primer set - knitting command will complain about
    # directory already existing
    mkdir $outDir . $primerFName . "_" . $primerRName . "/";

    # write primer input file to individual directory
    open my $individualPrimerFileFH, ">", $outDir . $primerFName . "_" . $primerRName . '/primerInfile.txt', or die "cannot create individual primer input file";
    print $individualPrimerFileFH $input, "\n";

    close $individualPrimerFileFH;

    # Write out qsub shell script for the individual primer set
    open my $qsubShellFH, ">", $outDir . $primerFName . "_" . $primerRName . '/qsub.sh', or die "cannot create qsub shell file\n";

    print $qsubShellFH $qsubHeader, "\n";

    # print out portion of header specific to primer
    print $qsubShellFH 
	'#$ -o ' . $outDir . $primerFName . "_" . $primerRName . '/qsubStdOut.txt' . "\n" .
	'#$ -e ' . $outDir . $primerFName . "_" . $primerRName . '/qsubStdErr.txt' . "\n" .
	'#$ -N PE_' .  $primerFName . "_" . $primerRName . "\n\n";

    # Specify R version to use -- Error here: use: "command not found" -- use -V in the qsub command and make sure to use r-3.5.1 before
    #print $qsubShellFH "source /usr/local/packages/sge-root/igs/common/settings.sh\n\nuse R-3.5.1\n\n";

    # R command
    my $rCommand = 'R -e \'rmarkdown::render("primerEvalPipeline.Rmd", ' . 
	'output_dir = "' . $outDir . $primerFName . "_" . $primerRName . '/", ' .
	'params = list(' .
	'outDir = "' . $outDir . $primerFName . "_" . $primerRName . '/", ' .
	'primerFile = "' . $outDir . $primerFName . "_" . $primerRName . '/primerInfile.txt", ' .
	'blastLoc = "' . $blastLoc . '", ' .
	'blastDb = "' . $blastDb . '", '  .
	'taxonomyDb = "' . $taxonomyDb . '", ' . 
	'threads = ' . $proc . '))\'';

    print $qsubShellFH $rCommand, "\n";

    close $qsubShellFH;

    # Run qsub command here
    my $qsubCommand = "qsub " . $outDir . $primerFName . "_" . $primerRName . '/qsub.sh';
    `$qsubCommand`;
}



##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    batchQsub.pl - submit qsub shell scripts for each primer in a multi-primer file
    
Usage:

    perl batchQsub.pl then put your arguments and stuff here


=head OPTIONS

Options:

    --verbose
    --help
    --input
    --outDir

=cut
