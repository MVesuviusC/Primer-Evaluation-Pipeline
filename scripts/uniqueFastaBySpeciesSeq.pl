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
my $taxa;
my $fasta;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
	    "taxonomy=s"        => \$taxa,
	    "fasta=s"           => \$fasta

      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $taxaHeader;
my @taxaCols;
my $giCol;
my $speciesCol;
my %taxaHash;
my %seqHash;

##############################
# Code
##############################


##############################
### Pull in taxa file and make hash -> $hash{gi} = species
### More stuff

open TAXAFILE, "$taxa" or die "Could not open taxonomy input\nWell, crap\n";
$taxaHeader = <TAXAFILE>;
chomp $taxaHeader;
@taxaCols = split "\t", $taxaHeader;

for(my $i = 0; $i < scalar(@taxaCols); $i++) {
    $giCol = $i if $taxaCols[$i] eq "gi";
    $speciesCol = $i if $taxaCols[$i] eq "species";
}

while (my $input = <TAXAFILE>){
    chomp $input;
    my @cols = split "\t", $input;
    if($cols[$speciesCol] !~ /sp\.|isol\.|isolate|cf\./) {
	$taxaHash{$cols[$giCol]} = $cols[$speciesCol];
    } 
}



##############################
### Pull in fasta file and keep one entry for each species:sequence combo
### More stuff
local $/ = "\n>";
open FASTAFILE, "$fasta" or die "Could not open fasta input\nWell, crap\n";
while (my $input = <FASTAFILE>){
    chomp $input;
    $input =~ s/^>//;
    my ($header, @rawSequence) = split "\n", $input;
    my $seq = join("", @rawSequence);
    $seq =~ s/[\t\s]//g;
    $seq = uc($seq);
    if(exists($taxaHash{$header})) {
	my $species = $taxaHash{$header};
	if(!exists($seqHash{$species . $seq})) {
	    print ">", $header, "\n", $seq, "\n";
	    $seqHash{$species . $seq} = 1;
	}
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
