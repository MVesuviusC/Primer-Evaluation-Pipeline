#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 10/28/19
# Last modified: 10/28/19
# Title: filterFastaNs.pl
# Purpose: Get rid of fasta entries with Ns in seqs
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $fastaIn;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
	    "fasta=s"           => \$fastaIn            
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################


##############################
# Code
##############################


##############################
### Stuff
### More stuff

$/ = "\n>";
open my $fastaFH, "$fastaIn" or die "Could not open fasta input\nWell, crap\n";

while (my $input = <$fastaFH>){
    chomp $input;
    my ($header, @sequences) = split "\n", $input;
    my $seq = join("", @sequences);

    $header =~ s/>//; # get rid of ">" in first entry

    if($seq !~ /N/i) {
	print ">", $header, "\n", $seq, "\n";
    }
}

$/ = "\n";



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
