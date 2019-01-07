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
my $input;
my $maxAmpLen;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "input=s"           => \$input,
            "maxAmpLen=i"       => \$maxAmpLen
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $lastSgi = "";
my $lastLine = "";

##############################
# Code
##############################


##############################
### Stuff
### More stuff

while (my $input = <>){
    chomp $input;
    my ($qseqid, $sgi, $qlen, $qstart, $qend, $sstart, $send, $slen, 
        $sstrand, $sscinames, $scomnames, $qseq, $sseq) 
            = split "\t", $input;
            
    if($sgi eq $lastSgi) {
        my ($lastQseqid, $lastSgi, $lastQlen, $lastQstart, $lastQend, 
            $lastSstart, $lastSend, $lastSlen, $lastSstrand, $lastSscinames, 
            $lastScomnames, $lastQseq, $lastSseq) 
                = split "\t", $lastLine;
                
        my @hitPositions = sort { $a <=> $b } ($sstart, $send, $lastSstart, $lastSend);
        
        if($sstrand eq $lastSstrand && ($hitPositions[3] - $hitPositions[0] < $maxAmpLen)) {
            # count primer mismatches
            my $lastMismatchCount = mismatchCounter($lastQseq, $lastSseq, $lastQstart);
            my $mismatchCount = mismatchCounter($qseq, $sseq, $qstart);
            
            print join("\t", $qseqid, $sgi, $sscinames, $scomnames, $qseq, $sseq, 
                $lastQseq, $lastSseq, $hitPositions[3] - $hitPositions[0], 
                $mismatchCount, $lastMismatchCount), "\n";
        }
    }
    $lastLine = $input;
    $lastSgi = $sgi;
}


sub mismatchCounter {
    my $seq1 = $_[0];
    my $seq2 = $_[1];
    my $seqPos = $_[2];
    
    my $mismatches = 0;
    
    for(my $i = 0; $i < length($seq1); $i++) {
        my $qBase = substr($seq1, $i, 1);
        my $sBase = substr($seq2, $i, 1);
        if($qBase ne $sBase) {
           $mismatches++;
        }
    }
    
    if($seqPos < 20) {
        $mismatches = "for:" . $mismatches;
    } else {
        $mismatches = "rev:" . $mismatches;
    }

    return $mismatches;
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
