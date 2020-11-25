#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 9-20-17
# Last modified: 11-12-20 
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
my $primerF;
my $primerR;
my $alignVar = 10;
my $bannedWords = "";

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "blastIn=s"         => \$blastIn,
            "taxaIn=s"          => \$taxaIn,
            "primerF=s"         => \$primerF,
            "primerrs"          => \$primerR,
            "alignVar=i"        => \$alignVar,
            "bannedWords=s"     => \$bannedWords
) or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);

##############################
# Global variables
##############################
my $forLength;
my $revLength;
my %taxaHash;
my $taxaCount = 0;
my $blastCount = 0;
my $notFoundTaxaCount = 0;
my $notFoundTaxaList;
my %resultsHash;

##############################
# Code
##############################

##############################
### Check that proper input was provided

if(!defined($primerF) || !defined($blastIn) || !defined($taxaIn)) {
    print STDERR "\n\nPrimer sequences and other input must be provided.\n";
    print STDERR "primerF: ", $primerF, "\nblastIn: ", $blastIn, "\ntaxaIn: ", $taxaIn, "\n\n";
    pod2usage(1);
    die;
}

##############################
### Read in primer sequences and store sequence lengths


$forLength = length($primerF);
$revLength = length($primerR);

##############################
### Pull in taxa info and make hash of form:
### $taxaHash{gi}{taxaLevel or taxId} = taxon

open my $taxaFH, "$taxaIn" or die "Could not open taxonomy input\n";
my @header = split "\t", <$taxaFH>;
chomp @header;
while(my $input = <$taxaFH>) {
    chomp $input;
    if($verbose) {
        $taxaCount++;
        print STDERR "Taxa entries processed: ", $taxaCount, "                      \r";
    }

    my ($taxid, $species, $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $tax_name) = split "\t", $input;

    # all entries have tax_name, but not all have species
    if($species eq "NA") {
        $species = $tax_name;
    }
    #construct taxa hash
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
### Go through blast output
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
        my ($query, $sTaxid, $score, $alignLen, $qStart, $qEnd, $qLen, $sStart, $sEnd, $sLen, $sAcc) = split "\t", $input;

        my $alignPercent = 100 * ($alignLen / $qLen);
    
      # check if alignment length is above a provided cutoff to avoid sequences with 30bp aligned out of a 300bp fragment
        # alignVar is a percent
        # there is a problem here. The output does not indicate which side each primer belongs on. For now I am assuming
        # the forward is on the front, but I may just change $forLength and $revLength to a numeric value option
        if($alignPercent >= (100 - $alignVar) && $alignPercent <= (100 + $alignVar)) {
          if($sStart >= ($forLength + $qStart - 1) && $sLen > ($sEnd + $revLength + ($qLen - $qEnd))) {
              # Put results into hash to remove duplicates
              if(exists($taxaHash{$sTaxid})) {
                  if($taxaHash{$sTaxid}{species} !~ /($bannedWords)/ || $bannedWords eq "") {
                    $resultsHash{$taxaHash{$sTaxid}{superkingdom} . "\t" . 
                            $taxaHash{$sTaxid}{kingdom} . "\t" . 
                            $taxaHash{$sTaxid}{phylum} . "\t" . 
                            $taxaHash{$sTaxid}{class} . "\t" . 
                            $taxaHash{$sTaxid}{order} . "\t" . 
                            $taxaHash{$sTaxid}{family} . "\t" . 
                            $taxaHash{$sTaxid}{genus} . "\t" . 
                            $taxaHash{$sTaxid}{species}
                        } = 1;
                  }
              } else {
                  $notFoundTaxaList .= $sTaxid . ", "; 
              print STDERR "taxid ", $sTaxid, " taxonomic information not found!\n";
                  $notFoundTaxaCount++;
              }
          }
      }
    }
}

if($notFoundTaxaCount > 0) {
    print STDERR "\n\n######\nA total of ", $notFoundTaxaCount, " entries did not have matching taxonomic information in the taxa input file.\n######\n\n";
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

=item B<--primerFile>

Sequences of the primers used.

=back

=cut
