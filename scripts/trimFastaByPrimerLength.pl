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
my $fastaInput;
my $primers;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "fastaInput=s"	=> \$fastaInput,
            "primers=s"	        => \$primers
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %primerHash;

##############################
# Code
##############################

##############################
### Get sequences to trim and deconvolute ambiguous bases

open PRIMERS, $primers or die "Could not open file of primer sequences to trim\n";

while (my $input = <PRIMERS>){
    chomp $input;
    my ($name1, $forward, $name2, $reverse) = split "\t", $input;
    $primerHash{$name1 . "_" . $name2}{forward} = length($forward);
    $primerHash{$name1 . "_" . $name2}{reverse} = (-1 * length($reverse));
}

##############################
### 

$/ = "\n>";
open SEQS, $fastaInput or die "Could not open fasta input\n";
while (my $input = <SEQS>){
    chomp $input;
    my $filename = $fastaInput;
    $filename =~ s/_sequences.fasta//;
    $filename =~ s/.+\///;
    my ($header, @seqArray) = split "\n", $input;
    $header =~ s/^>//; # replace ">" in first fasta entry
    my $sequence = join("", @seqArray);
    $sequence =~ s/[\s\t-]//g;
    $sequence = uc($sequence);
    if($sequence ne "NA") {
	$sequence = substr($sequence, $primerHash{$filename}{forward});
	$sequence = substr($sequence, 0, $primerHash{$filename}{reverse});

	print ">", $header, "\n", $sequence, "\n";
    }
}
$/ = "\n";
close SEQS;

###
##############################

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
