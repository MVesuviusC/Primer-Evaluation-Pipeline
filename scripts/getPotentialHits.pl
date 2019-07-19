#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 9-20-17
# Last modified: 9-20-17 
# Title: getPotentialHits.pl
# Purpose: Parse BLAST hits
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $blastIn;
my $taxaIn;
my $primerF = "";
my $primerR = "";
my $alignVar = 0;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "blastIn=s"		=> \$blastIn,
            "taxaIn=s"		=> \$taxaIn,
	    "primerF=s"         => \$primerF,
	    "primerR=s"         => \$primerR,
	    "alignVar=i"        => \$alignVar
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %taxaHash;
my $taxaCount = 0;
my $blastCount = 0;
my $notFoundTaxaCount = 0;
my $notFoundTaxaList;
my $forLength = length($primerF);
my $revLength = length($primerR);
my %resultsHash;

##############################
# Code
##############################

##############################
### Check that proper input was provided

if($primerF eq "" || $primerR eq "" || !defined($blastIn) || !defined($taxaIn)) {
    print STDERR "\n\nPrimer sequences and other input must be provided.\n";
    print STDERR "primerF: ", $primerF, "\nprimerR: ", $primerR, "\nblastIn: ", $blastIn, "\ntaxaIn: ", $taxaIn, "\n\n";
    pod2usage(1);
    die;
}


##############################
### Pull in taxa info and make hash of form:
### $taxaHash{gi}{taxaLevel or taxId} = taxon

open TAXA, "$taxaIn" or die "Could not open taxonomy input\n";
my @header = split "\t", <TAXA>;
chomp @header;
while(my $input = <TAXA>) {
    chomp $input;
    if($verbose) {
	$taxaCount++;
	print STDERR "Taxa entries processed: ", $taxaCount, "                      \r";
    }
    $input =~ s/\n//g;
    my @columns = split "\t", $input;
    for(my $i = 1; $i < scalar(@columns); $i++) {
	if($header[$i] eq "species") {
	    $columns[$i] =~ s/.+\ssp\..*/NA/;
	    $columns[$i] =~ s/.*uncultured.*/NA/;
	    $columns[$i] =~ s/.*unidentified.*/NA/;
	    $columns[$i] =~ s/.*symbiont.*/NA/;
	    $columns[$i] =~ s/.+\scf\.\s.+/NA/;
	    $columns[$i] =~ s/.+isolate\s.+/NA/;
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
 
# want: query subjectGi subjectOrder subjectFamily subjectGenus subjectSpecies

#print "query\tsubjectGi\tqueryKingdom\tqueryPhylum\tqueryClass\tqueryOrder\tqueryFamily\tquerySpecies\tsubjectOrder\tsubjectFamily\tsubjectGenus\tsubjectSpecies\n";

while (my $input = <BLAST>) {
    chomp $input;
    
    if($input !~ /^\#/) {
	if($verbose) {
	    $blastCount++;
	    print STDERR "Blast entries processed: ", $blastCount, "                     \r";
	}
	## Need: qlen, slen to figure out if the whole thing aligned or if just part and if the subject has additional sequence that the
	## primers could have hit
	my ($query, $subject, $evalue, $alignLen, $qLen, $sStart, $sEnd, $sLen) = split "\t", $input;

	
	my $alignPercent = 100 * ($alignLen / $qLen);

	if($alignPercent >= (100 - $alignVar) && $alignPercent <= (100 + $alignVar)) {      #$qLen == $alignLen) {
	    if(($sStart - $forLength) >= 0 && ($sLen - $sEnd) >= $revLength ) {
		# Put results into hash to remove duplicates
		if(exists($taxaHash{$subject})) {
		    $resultsHash{$taxaHash{$subject}{superkingdom} . "\t" . 
				     $taxaHash{$subject}{kingdom} . "\t" . 
				     $taxaHash{$subject}{phylum} . "\t" . 
				     $taxaHash{$subject}{class} . "\t" . 
				     $taxaHash{$subject}{order} . "\t" . 
				     $taxaHash{$subject}{family} . "\t" . 
				     $taxaHash{$subject}{genus} . "\t" . 
				     $taxaHash{$subject}{species}
			     } = 1;
		} else {
		    $notFoundTaxaList .= $subject . ", "; 
                    #print STDERR "gi# ", $subject, " taxonomic information not found!\n";
		    $notFoundTaxaCount++;
		}
	    }
	}
    }
}

if($notFoundTaxaCount > 0) {
    print STDERR "\n\n######\nA total of ", $notFoundTaxaCount, " gis did not have matching taxonomic information in the taxa input file.\n######\n\n";
    #print STDERR "Gi list:", $notFoundTaxaList, "\n";
}

##############################
### Print out results

if($verbose) {
    print STDERR "Printing output\n";
}

print "superkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";

for my $output (keys %resultsHash) {
    print $output, "\n";
}

##############################
# POD
##############################

    
=head SYNOPSIS

Summary:    
    
    getPotentialHits.pl - Using blast output, taxa information and primer sequences, print out any blast hits that could have included the primer sequences
    
Usage:

    perl getPotentialHits.pl --blastIn <blast.outfmt7> --taxaIn <taxa.txt> --primerF <forward primer> --primerR <reverse primer> [options] 


=head OPTIONS

Options:

=over 4

=item B<--verbose>

Print out updates as the program runs.

=item B<--help>

Print out some "helpful" information

=item B<--taxaIn>

Taxonomic information on BLAST hits. Should be the gi number followed by columns of the taxonoic information. Needs to include a header line. 

=item B<--blastIn>

Blast input. Should be output from BLAST with -outfmt "7 qseqid sgi length qlen sstart send slen" 

=item B<--primerF>

Sequence of the forward primer used.

=item B<--primerR>

Sequence of the reverse primer used.

=back

=cut
