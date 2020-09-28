#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 
# Last modified: 
# Title: countNumInRank.pl
# Purpose: 
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $input;
my $outDir = '.';

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
      	    "input=s"           => \$input,
      	    "outDir=s"          => \$outDir
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);

##############################
# Global variables
##############################
my %taxaStorageHash;

##############################
# Code
##############################
$outDir .= '/';

##############################
### Store all data in a hash to keep unique entries

# need to kick out any blast hit containing "banned words" 
#     (sp., cf., isolate, uncultured, symbiont, unidentified)
# should pull these in from a file

open my $inputFH, "$input" or die "Could not open input\nWell, crap\n";

my $header = <$inputFH>;
chomp $header;
$header =~ s/subject//g; # cut the label to just taxa name
my @labelArray = split "\t", $header;

while (my $input = <$inputFH>){
    chomp $input;
    my @inputArray = split "\t", $input;
    my $query = $inputArray[0];
	
	  # order through species
	  for(my $level = 7; $level <=  10; $level++) {
	    if($inputArray[$level] ne "NA") {
    		# Hash organization: $taxaStorageHash{Taxa level}{query ID}{Taxa name}
    		# This keeps one copy of each taxa name for each query for each taxa level so I can count later
		    $taxaStorageHash{$labelArray[$level]}{$query}{$inputArray[$level]} = 1;
	    }	
	  }
}

##############################
### Count number of taxa names per level per query and print

open my $longOutputFH, ">", $outDir . "topHitSummary.txt";
open my $averageOutputFH, ">", $outDir . "topHitMeans.txt";

print $longOutputFH "TaxaLevel\tquery\tUniqueTaxaHit\n";
print $averageOutputFH "TaxaLevel\tMeanNumTaxaHit\n";

for my $taxaLevel (keys %taxaStorageHash) {
    my $taxaCountSum = 0;
    my $queryNum = 0;
    for my $queryID (keys %{ $taxaStorageHash{$taxaLevel} }) {
    	my @keyArray = keys %{ $taxaStorageHash{$taxaLevel}{$queryID} };
    	my $numTaxa = scalar(@keyArray);
    
    	print $longOutputFH $taxaLevel, "\t", $queryID, "\t", $numTaxa, "\n";
    	
    	if($numTaxa != 0) { # skip the queries that hit nothing useful
    	    $queryNum++;
    	    $taxaCountSum += $numTaxa;
    	}
    }
    if($queryNum != 0) {
    	print $averageOutputFH $taxaLevel, "\t", $taxaCountSum / $queryNum, "\n";
    } else {
    	print $averageOutputFH $taxaLevel, "\tND\n";
    }
}





##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    countNumInRank.pl - counts the number of taxonomic groups hit by each query
    
Usage:

    perl countNumInRank.pl [options] --input inFile.txt > outFile.txt


=head OPTIONS

Options:

    --verbose
    --help
    --input
    --outDir

=cut
