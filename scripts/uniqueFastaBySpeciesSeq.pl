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
my $fasta;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
	    "fasta=s"           => \$fasta
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %seqHash;

##############################
# Code
##############################


##############################
### Pull in fasta file and keep one entry for each species:sequence combo

local $/ = "\n>";
open FASTAFILE, "$fasta" or die "Could not open fasta input\nWell, crap\n";
while (my $input = <FASTAFILE>){
    chomp $input;
    $input =~ s/^>//;
    my ($header, @rawSequence) = split "\n", $input;
    my $seq = join("", @rawSequence);
    $seq =~ s/[\t\s]//g;
    $seq = uc($seq);

    my $species = $header;
    $species =~ s/:sk-.+//; # cut off everything but the species
    
    if(!exists($seqHash{$species . $seq})) {
	print ">", $header, "\n", $seq, "\n";
	$seqHash{$species . $seq} = 1;	
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
