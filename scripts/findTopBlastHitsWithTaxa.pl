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
my %taxaHash;
my $lastQuery = "";
my $minEVal = 0;

##############################
# Code
##############################

##############################
### Pull in taxa info and make hash of form:
### $taxaHash{gi}{taxaLevel or taxId} = taxon

open TAXA, "$taxaIn" or die "Could not open taxonomy input\n";
my @header = split "\t", <TAXA>;
chomp @header;
while(my $input = <TAXA>) {
    chomp $input;
    $input =~ s/\n//g;
    my @columns = split "\t", $input;
    for(my $i = 1; $i < scalar(@columns); $i++) {
	if($header[$i] eq "species") {
	    $columns[$i] =~ s/.+\ssp\..*/NA/;
	    $columns[$i] =~ s/.+\scf\.\s.+/NA/;
	    $columns[$i] =~ s/.*isolate\s.*/NA/;
            $columns[$i] =~ s/.*uncultured.*/NA/;
            $columns[$i] =~ s/.*unidentified.*/NA/;
            $columns[$i] =~ s/.*symbiont.*/NA/;

	    $columns[$i] =~ s/\s/_/;
	    $columns[$i] =~ s/\s.+//;
	}
	$taxaHash{$columns[0]}{$header[$i]} = $columns[$i]; 
    }
}

##############################
### Go through blast output and print out top hits with taxa
### 

$blastIn =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
open BLAST, "$blastIn" or die "Could not open BLAST input\n";

print "query\tsubjectGi\tqueryKingdom\tqueryPhylum\tqueryClass\tqueryOrder\tqueryFamily\tquerySpecies\tsubjectOrder\tsubjectFamily\tsubjectGenus\tsubjectSpecies\n";

while (my $input = <BLAST>) {
    chomp $input;
    if($input !~ /#/) {
	my ($query, $subject, $evalue, $alignLen, $qLen, $sStart, $sEnd, $sLen) = split "\t", $input;
	# if new query, print out min eval and gis and reset variables
	if($lastQuery ne $query) {
	    $minEVal = 1;
	}
	# compare minEval to evalue and keep if equal or new 
	if($evalue <= $minEVal) {
	    my $sGi = $subject;
	    $sGi =~ s/^gi.//;
	    $sGi =~ s/\|.+//;
	    if(exists($taxaHash{$sGi}) && exists($taxaHash{$query})) {
		my $genus = $taxaHash{$sGi}{species};
		$genus =~ s/_.+//;
		#print STDERR $query, "\n";
		print join("\t", 
			   $query, 
			   $sGi,
			   $taxaHash{$query}{kingdom},
			   $taxaHash{$query}{phylum},
			   $taxaHash{$query}{class},
			   $taxaHash{$query}{order},
			   $taxaHash{$query}{family},
			   $taxaHash{$query}{species},
			   $taxaHash{$sGi}{order}, 
			   $taxaHash{$sGi}{family}, 
			   $genus, 
			   $taxaHash{$sGi}{species},
		    ), "\n";
	    }
	    $minEVal = $evalue;
	}
	$lastQuery = $query;
    }
}


##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    xxxxxx.pl - generates a consensus for a specified gene in a specified taxa
    
Usage:

    perl xxxxxx.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
