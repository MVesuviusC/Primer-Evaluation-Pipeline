#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 
# Last modified: 
# Title: .pl
# Purpose: 
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $taxaIn;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "taxaIn=s"		=> \$taxaIn
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %taxaHash;

##############################
# Code
##############################

##############################
### Pull in taxa info and make hash

open my $taxaFh, "$taxaIn" or die "Could not open taxonomy input\n";
my @header = split "\t", <$taxaFh>;
chomp @header;
while(my $input = <$taxaFh>) {
    chomp $input;
    my ($taxid, $species, $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $tax_name) = split "\t", $input;

    if($species eq "NA") {
	$species = $tax_name;
    }
    my $fullTaxonomy = join("\t", $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $species);

    $taxaHash{$fullTaxonomy} = 1;
}

##############################
### Go through hash and print out list of taxa

for my $line (sort keys %taxaHash) {
    print $line, "\n";
}



##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    xxxxxx.pl - 
    
Usage:

    perl xxxxxx.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
