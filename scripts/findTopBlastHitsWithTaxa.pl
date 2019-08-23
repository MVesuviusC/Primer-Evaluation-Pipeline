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
my $blastIn;
my $taxaIn;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "blastIn=s"		=> \$blastIn,
            "taxaIn=s"		=> \$taxaIn
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $fullQueryName;
my %taxaHash;
my $maxScore = 0;

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

    $taxaHash{$taxid}{species} = $species;
    $taxaHash{$taxid}{superkingdom} = $superkingdom;
    $taxaHash{$taxid}{kingdom} = $kingdom;
    $taxaHash{$taxid}{phylum} = $phylum;
    $taxaHash{$taxid}{class} = $class;
    $taxaHash{$taxid}{order} = $order;
    $taxaHash{$taxid}{family} = $family;
    $taxaHash{$taxid}{genus} = $genus;
}

##############################
### Go through blast output and print out top hits with taxa
### 

$blastIn =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
open my $blastFh, "$blastIn" or die "Could not open BLAST input\n";

print "query\tsubjectAccession\tsubjectTaxid\tsubjectSuperkingdom\tsubjectKingdom\tsubjectPhylum\tsubjectClass\tsubjectOrder\tsubjectFamily\tsubjectGenus\tsubjectSpecies\n";

while (my $input = <$blastFh>) {
    chomp $input;
    if($input !~ /^\#/) {
	my ($query, $sTaxid, $score, $alignLen, $qLen, $sStart, $sEnd, $sLen, $sAcc) = split "\t", $input;
	# compare minEval to evalue and keep if equal or new 
	if($score >= $maxScore) {
	    if(exists($taxaHash{$sTaxid})) {
		# Don't print any hit that has uncertain taxonomic assignment
		if($taxaHash{$sTaxid}{species} !~ /sp\.|cf\.|isolate|uncultured|symbiont|unidentified|NA/i) {
		    print join("\t", 
			       $fullQueryName, 
			       $sAcc,
			       $sTaxid,
			       $taxaHash{$sTaxid}{superkingdom},
			       $taxaHash{$sTaxid}{kingdom},
			       $taxaHash{$sTaxid}{phylum},
			       $taxaHash{$sTaxid}{class},
			       $taxaHash{$sTaxid}{order}, 
			       $taxaHash{$sTaxid}{family}, 
			       $taxaHash{$sTaxid}{genus}, 
			       $taxaHash{$sTaxid}{species},
			       ), "\n";
		}
	    }
	    $maxScore = $score;
	}
    } else {
	if($input =~ /Query/) {
	    $fullQueryName = $input;
	    $fullQueryName =~ s/\# Query: //;
	}
	# if new query, reset maxScore
	$maxScore = 0;
    }
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
