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
       A => ["A"],
       T => ["T"],
       G => ["G"],
       C => ["C"],
       W => ["A","T"],
       S => ["C","G"],
       M => ["A","C"],
       K => ["G","T"],
       R => ["A","G"],
       Y => ["C","T"],
       B => ["C","G","T"],
       D => ["A","G","T"],
       H => ["A","C","T"],
       V => ["A","C","G"],
       N => ["A","C","G","T"]
       );

my %degenerateRegexHash = ( #hash of arrays - degenerate bases with matching bases
       A => "A",
       T => "T",
       G => "G",
       C => "C",
       W => "AT",
       S => "CG",
       M => "AC",
       K => "GT",
       R => "AG",
       Y => "CT",
       B => "CGT",
       D => "AGT",
       H => "ACT",
       V => "ACG",
       N => "ACGT"
       );

my %primerHash; # primerHash{primerName} = (forward, reverse)
my $nCount = 20;
my $lastSgi = "";
my $lastLine = "";
my $blastEntryCounter = 1;

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
    
    #### Need to build in check to be sure format is right and force uppercase
    my ($primerName, $primerF, $primerR) = split "\t", $input;
    
    # Make up hash to use for parsing blast output
    if(exists($primerHash{$primerName})) {
        die "Primer sets must have unique names\n";
    } else {
        @{ $primerHash{$primerName} } = ($primerF, revComp($primerR));
    }
    
    # deconvolute ambiguious bases for both primers
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
            print SEQSOUT ">", $primerName, "\n", $seq1, "N" x $nCount, revComp($seq2), "\n";
            $i++;
        }
    }
}

close SEQSOUT; 

if($verbose) {
    print STDERR "starting BLAST\n";
}

### need to get taxid from blast results
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
    "-outfmt \"6 qseqid sgi qlen qstart qend sstart send slen sstrand sscinames scomnames qseq sseq staxids\" " .
    "-ungapped " .
    "-sum_stats true";


open BLASTRESULTS, "-|", $blastCmd or die "Blast query failed\n";

#open BLASTRESULTS, "testBlast.txt" or die "crap\n";

while(my $input = <BLASTRESULTS>) {
    chomp $input;
    my @parsed = parseBlast($input);
    if(scalar(@parsed) > 0) {
        ########### Need to get rid of redundant hits due to ambiguous bases
        my $sgi = 
        print join("\t", @parsed), "\n";
    }
    if($verbose) {
        print STDERR "Parsing BLAST entry $blastEntryCounter                  \r";
        $blastEntryCounter++;
    }
}

system("rm $tempName" . ".txt ");

if($verbose) {
    print STDERR "\nDone!\n";
}


######################
### Subfunctions

sub addDegeneratePrimer {
    my $primer = shift;
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
    return uniq(@primerArray);
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

sub revComp {
    my $seq = shift;
    $seq =~ tr/ACGTacgtYRWSKMDVHBXN/TGCAtgcaRYWSMKHBDVXN/;
    $seq = reverse ($seq);
    return $seq;
}

sub uniq {
    my @array = @_;
    my %uniqHash;
    for my $element (@array) {
        $uniqHash{$element} = 1;
    }
    return keys %uniqHash;
}

sub parseBlast {
    my $blastData = shift;
    my ($qseqid, $sgi, $qlen, $qstart, $qend, $sstart, $send, $slen,
        $sstrand, $sscinames, $scomnames, $qseq, $sseq, $staxids)
            = split "\t", $blastData;
    my @output;
    
    if($sgi eq $lastSgi) {
        my ($lastQseqid, $lastSgi, $lastQlen, $lastQstart, $lastQend,
            $lastSstart, $lastSend, $lastSlen, $lastSstrand, $lastSscinames,
            $lastScomnames, $lastQseq, $lastSseq, $lastStaxids)
                = split "\t", $lastLine;

        my @hitPositions = sort { $a <=> $b } ($sstart, $send, $lastSstart, $lastSend);

        if($sstrand eq $lastSstrand && ($hitPositions[3] - $hitPositions[0] < $maxAmpLen)) {
            # count primer mismatches
            my $lastMismatchCount = mismatchCounter($lastLine); #$lastQseq, $lastSseq, $lastQstart, $lastQend, $lastQseqid, $lastSgi, $lastSstart);
            my $mismatchCount = mismatchCounter($blastData); #$qseq, $sseq, $qstart, $qend, $qseqid, $sgi);
            @output = (
                        $qseqid,                                # query ID
                        $sgi,                                   # hit gi
                        $staxids,                               # taxid of hit - can be multiple
                        $sscinames,                             # scientific name
                        $scomnames,                             # common name
                        $hitPositions[0],                       # start of amplicon location
                        $hitPositions[3],                       # end amplicon location
                        $hitPositions[3] - $hitPositions[0],    # amplicon length
                        $mismatchCount,                         # number of mismatches
                        $lastMismatchCount,                     # number of mismatches other end
                        $qseq,                                  # query seq
                        $sseq,                                  # hit seq
                        $lastQseq,                              # other end query seq
                        $lastSseq                               # other end hit seq
                      );
        }
    }
    $lastLine = $blastData;
    $lastSgi = $sgi;
    return @output;
}


sub mismatchCounter {
    my $blastData = shift;
    my ($qseqid, $sgi, $qlen, $qstart, $qend, $sstart, $send, $slen,
        $sstrand, $sscinames, $scomnames, $qseq, $sseq, $staxids)
            = split "\t", $blastData;

    my $mismatches = 0;
    #my $mismatchLocs; Will have to figure out how best to implement this later.... 
    
    # Count the mismatches in the aligned portion of the primer
    for(my $i = 0; $i < length($qseq); $i++) {
        my $qBase = substr($qseq, $i, 1);
        my $sBase = substr($sseq, $i, 1);
        if($qBase ne $sBase) {
           $mismatches++;
           #$mismatchLocs .= $sBase;
        } #else {
            #$mismatchLocs .= ".";
        #}
    }
    
    # Check if the full length of the primer is aligned and count additional mismatches if not
    if($qstart < 20) { # query is the forward primer
        my $primerShouldStart = 1;
        my $primerShouldEnd = length(@{ $primerHash{$qseqid} }[0]);
        #print STDERR $primerShouldStart, "\t", $primerShouldEnd, "\t", $qstart, "\t", $qend, "\n", $blastData, "\n";
        if($qstart != $primerShouldStart || $qend != $primerShouldEnd) {
          #print STDERR "too short\n";
          my $primerSeq = @{ $primerHash{$qseqid} }[0];
          $mismatches += addUnalignedMismatches($blastData, $primerSeq, $primerShouldStart, $primerShouldEnd); #, $mismatchLocs);
        }
        $mismatches = "for:" . $mismatches;
    } else { # query is the reverse primer
        my $primerShouldStart = length(@{ $primerHash{$qseqid} }[0]) + $nCount + 1;
        my $primerShouldEnd = $primerShouldStart + length(@{ $primerHash{$qseqid} }[1]) - 1;

        if($qstart != $primerShouldStart || $qend != $primerShouldEnd) {
          my $primerSeq = @{ $primerHash{$qseqid} }[1];
          $mismatches += addUnalignedMismatches($blastData, $primerSeq, $primerShouldStart, $primerShouldEnd); #, $mismatchLocs);
        }
        $mismatches = "rev:" . $mismatches;
    }

    return $mismatches;
}

sub addUnalignedMismatches {
    my $blastData = $_[0];
    my $primerSeq = $_[1];
    my $primerShouldStart = $_[2];
    my $primerShouldEnd = $_[3];
    #my $mismatchLocString = $_[4];
    
    my ($qseqid, $sgi, $qlen, $qstart, $qend, $sstart, $send, $slen,
        $sstrand, $sscinames, $scomnames, $qseq, $sseq, $staxids)
            = split "\t", $blastData;

    my $mismatchCount = 0;
    
    #Check 5' end of primer
    if($qstart == ($primerShouldStart + 1)) { # only first base is a mismatch
        $mismatchCount++;
        #$mismatchLocString = "." . $mismatchLocString;
    } elsif($qstart > ($primerShouldStart + 1)) { # more than one mismatch at 5' end of forward primer
        $mismatchCount++; # increment because we know the base at $qstart - 1 is wrong
        #loop from 1 to $qStart - 1, get base with blastdbcmd and compare with primer seq
        for(my $i = 0; $i < ($qstart - $primerShouldStart - 1); $i++) { # don't want to compare base at $qstart - 1; we know it's wrong
            my $primerBase = substr($primerSeq, $i, 1); # -1 due to 0 start
            my $seqBaseNumberToGet = $sstart - ($qstart - $primerShouldStart) + $i;
            my $getSeqCmd = "blastdbcmd -outfmt \"%s\" -db " . $blastDb . " -entry " . $sgi . " -range " . $seqBaseNumberToGet;
            my $seqBase = `$getSeqCmd`;
            chomp $seqBase;
            $seqBase = uc($seqBase);
            if($seqBase !~ /[$degenerateRegexHash{$primerBase}]/) {
                $mismatchCount++;
            }
        }
    }
    
    #Check 3' end of primer
    if($qend == ($primerShouldEnd - 1)) { # only last base is a mismatch
        $mismatchCount++;
    } elsif($qend < ($primerShouldEnd - 1)) {
        $mismatchCount++; # increment because we know the base at $qend + 1 is wrong
        #loop from $qend + 1 to $primerShouldEnd, get base with blastdbcmd and compare with primer seq
        for(my $i = 1; $i < $primerShouldEnd - $qend; $i++) { # don't want to compare $qend + 1; we know it's wrong
            my $primerBase = substr($primerSeq, length($primerSeq) - ($primerShouldEnd - $qend) + $i, 1);
            my $seqBaseNumberToGet = $send + $i + 1;
            my $getSeqCmd = "blastdbcmd -outfmt \"%s\" -db " . $blastDb . " -entry " . $sgi . " -range " . $seqBaseNumberToGet;
            my $seqBase = `$getSeqCmd`;
            chomp $seqBase;
            $seqBase = uc($seqBase);
            if($seqBase !~ /[$degenerateRegexHash{$primerBase}]/) {
                $mismatchCount++;
            }
        }
    }
    return $mismatchCount;
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
