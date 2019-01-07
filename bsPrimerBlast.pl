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
my $primerInput;
my $blastDb;
my $processors;
my $tempName = "tempFile_" . int(rand(100000));
my $maxAmpLen = 2000;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "primerInput=s"     => \$primerInput,
            "blastDb=s"         => \$blastDb,
            "processors=i"      => \$processors,
            "tempName=s"        => \$tempName,
            "maxAmpLen=i"       => \$maxAmpLen
            
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %degeneratehash = ( #hash of arrays - degenerate bases with matching bases
                       W =>["A","T"],
                       S =>["C","G"],
                       M =>["A","C"],
                       K =>["G","T"],
                       R =>["A","G"],
                       Y =>["C","T"],
                       B =>["C","G","T"],
                       D =>["A","G","T"],
                       H =>["A","C","T"],
                       V =>["A","C","G"],
                       N =>["A","C","G","T"]
                       );


##############################
# Code
##############################


##############################
### Stuff
### More stuff
open SEQSOUT, ">", $tempName . ".txt";

open PRIMERINPUTFILE, "$primerInput" or die "Could not open primer input\nWell, crap\n";

while (my $input = <PRIMERINPUTFILE>) {
    chomp $input;
    my ($primerName, $primerF, $primerR) = split "\t", $input;
    # deconvolute ambiguious bases for both primers
    if($primerF =~ /[WSMKRYBDHVN]/ || $primerR =~ /[WSMKRYBDHVN]/) { #if the primer has any degenerate bases, deconvolute those and add the subsequent primers to the hash
        my @primerFArray = ($primerF);
        my @primerRArray = ($primerR);
        
        if($primerF =~ /[WSMKRYBDHVN]/) {
            @primerFArray = addDegeneratePrimer($primerF);
        }
        
        if($primerR =~ /[WSMKRYBDHVN]/) {
            @primerRArray = addDegeneratePrimer($primerR);
        }
        
        my $i = 1;
        for my $seq1 (@primerFArray) {
            for my $seq2 (@primerRArray) {
                print SEQSOUT ">", $primerName, "\n", $seq1, "N" x 20, revComp($seq2), "\n";
                $i++;
            }
        }
        
    } else {
        # write out to temp fasta file and make N-buffered sequence combinations
        print SEQSOUT ">", $primerName, "\n", $primerF, "N" x 20, revComp($primerR), "\n";
    }
}

close SEQSOUT; 

my $blastCmd =
    "blastn " .
    "-db " . $blastDb . " " .
    "-query " . $tempName . ".txt " .
    "-task blastn " .
    "-evalue 30000 " .
    "-word_size 7 " .
    "-max_hsps 100 " .
    "-max_target_seqs 50000 " .
    "-num_threads " . $processors . " " .
    "-reward 1 " .
    "-penalty -1 " .
    "-gapopen 2 " .
    "-gapextend 1 " .
    "-outfmt \"6 qseqid sgi qlen qstart qend sstart send slen sstrand sscinames scomnames qseq sseq\" " .
    "-ungapped " .
    "-sum_stats true " .
    "| " . "perl blastProcessing.pl --input - --maxAmpLen $maxAmpLen";

# maybe call a separate function to process the blast data as it comes in 
# instead of holding it all in memory... 
# need to get back primerName, gi, taxid, sscinames, scomnames, qseq, sseq, 
my $blastResults = `$blastCmd`;

print $blastResults;

system("rm $tempName" . ".txt ");

######################
### Subfunctions

sub addDegeneratePrimer {
    my $primer = $_[0];
    my @primerArray = ($primer); #make an array containing the degenerate primer
    my @tempArray = (); #make a temporary array to hold the new versions of the primers
    my $test = 1;
    while($test == 1){
        for(my $i = 0; $i < scalar(@primerArray); $i++) { # sort through primerArray
            if($primerArray[$i] =~ /[WSMKRYBDHVN]/) {
                push( @tempArray, getNewPrimerVersions($primerArray[$i]) );
            } else {
                push( @tempArray, $primerArray[$i] ); #add normal primer to tempArray
            }
        }
        @primerArray = @tempArray;
        @tempArray = ();
        if(join("", @primerArray) =~ /[WSMKRYBDHVN]/ == 0){
            $test = 0;
        }
    }
    return @primerArray;
}

sub getNewPrimerVersions {
    my $primer = shift; #primer sequence
    my @baseArray = split("", $primer); #split the primer up into individual bases in an array
    my @tempArray=();
    for(my $i = 0; $i < scalar(@baseArray); $i++){ #go through each base in the primer
        my $nucleotide = $baseArray[$i];
        if($nucleotide =~ /[WSMKRYBDHVN]/) { #if that base has a degenerate base
            for(@{$degeneratehash{$nucleotide}}) { #go through the possible replacements
                my @copyArray = @baseArray; # copy this array so we can modify it
                $copyArray[$i] = $_; #then switch out the base in the copy array with a possibility
                push(@tempArray, join("", @copyArray)); #and add the new decoded primer sequence to the end of the @primerArray
            }
        }
    }
    return @tempArray;
}

sub revComp{
    my $seq = shift;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    $seq = reverse ($seq);
    return $seq;
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
