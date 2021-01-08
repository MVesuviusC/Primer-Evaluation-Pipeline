#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 11/06/20
# Last modified: 11/06/20
# Title: prepTaxids.pl
# Purpose: Prep taxids for getTaxaLocal.pl
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $inFile;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "inFile=s"          => \$inFile
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %taxids;

##############################
# Code
##############################

$inFile =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;

##############################
### Stuff
### More stuff

open my $blastResultsInput, $inFile or die "Could not open input\nWell, crap\n";

while (my $input = <$blastResultsInput>){
    chomp $input;
    if($input !~ /\#/) {
        my @lineArray = split "\t", $input;
        # some "taxids" returned have multiple separated by ";"
        my @taxidArray = split ";", $lineArray[1];
        for(@taxidArray) {
          $taxids{$_} = 1;
        }
    }
}

for my $taxid (keys %taxids) {
    print $taxid, "\n";
}


##############################
# POD
##############################

#=pod

=head SYNOPSIS

Summary:

    filterFastaNs.pl - generates a consensus for a specified gene in a specified taxa

Usage:

    perl filterFastaNs.pl [options] --fasta infile.fasta


=head OPTIONS

Options:

    --verbose
    --help
    --fasta

=cut
