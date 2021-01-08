#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);

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
my $tempDir                 = ".";
my $maxPrimersPerFile       = 100;
my $maxPrimerVariantsTested = 2000;
my $forward;
my $reverse;
my $primerName = "primer";
my $blastVer   = "blastn";
my $blastDb;
my $processors;
my $tempName;
my $minAmpLen          = 1;
my $maxAmpLen          = 2000;
my $primerTipLen       = 5;
my $primerTipMismatch  = 2;
my $totalMismatchCount = 6;
my $debug;

# i = integer, s = string
GetOptions(
    "verbose"                   => \$verbose,
    "help"                      => \$help,
    "primerInput=s"             => \$primerInput,
    "tempDir=s"                 => \$tempDir,
    "maxPrimersPerFile=i"       => \$maxPrimersPerFile,
    "maxPrimerVariantsTested=s" => \$maxPrimerVariantsTested,
    "forward=s"                 => \$forward,
    "reverse=s"                 => \$reverse,
    "primerName=s"              => \$primerName,
    "blastVer=s"                => \$blastVer,
    "blastDb=s"                 => \$blastDb,
    "processors=i"              => \$processors,
    "tempName=s"                => \$tempName,
    "minAmpLen=i"               => \$minAmpLen,
    "maxAmpLen=i"               => \$maxAmpLen,
    "primerTipLen=i"            => \$primerTipLen,
    "primerTipMismatch=i"       => \$primerTipMismatch,
    "totalMismatch=i"           => \$totalMismatchCount,
    "debug"                     => \$debug
) or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);

##############################
# Global variables
##############################
my %degeneratehash = (    #hash of arrays - degenerate bases with matching bases
    A => ["A"],
    T => ["T"],
    G => ["G"],
    C => ["C"],
    W => [ "A", "T" ],
    S => [ "C", "G" ],
    M => [ "A", "C" ],
    K => [ "G", "T" ],
    R => [ "A", "G" ],
    Y => [ "C", "T" ],
    B => [ "C", "G", "T" ],
    D => [ "A", "G", "T" ],
    H => [ "A", "C", "T" ],
    V => [ "A", "C", "G" ],
    N => [ "A", "C", "G", "T" ]
);

my %degenerateRegexHash =
  (    #hash of arrays - degenerate bases with matching bases
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
    N => "ACGT",
    X => "ACGT"
  );

my %primerHash; # primerHash{primerName} = (forward, reverse) reverse is revComp
my $primerOutFileNum  = 1;
my $nCount            = 20;
my $lastSgi           = "";
my $lastLine          = "";
my $blastEntryCounter = 1;
my %resultsHash;
my $numberOfBlastcmdCalls = 0;

##############################
# Code
##############################

##############################
### Check inputs

# Check if -f/-r and --primerInput are provided and error out if TRUE
# if -f/-r, make sure --primerName is provided
# Check if Blast db provided

##############################
### Make temp file name if not defined

if ( !$tempName ) {
    my @numChars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
    $tempName = "tempFile_";
    $tempName .= $numChars[ rand @numChars ]
      for 1 .. 15;    # https://www.perlmonks.org/?node_id=233023
}

##############################
### Create file of primer(s) for blastn input

if ($forward) {       # primers provided directly
    addPrimersToHash( $forward, $reverse, $primerName );
}
else {                # list of primers provided through file
    open PRIMERINPUTFILE, "$primerInput"
      or die "Could not open primer input\nWell, crap\n";

    while ( my $input = <PRIMERINPUTFILE> ) {
        chomp $input;

        #### Need to build in check to be sure format is right and force uppercase
        #### Also check that just one primer set is passed
        my ( $primerNameF, $primerF, $primerNameR, $primerR ) = split "\t",
          $input;

        $forward = $primerF;
        $reverse = $primerR;

        $primerName = $primerNameF . "_" . $primerNameR;

        addPrimersToHash( $primerF, $primerR, $primerName );
    }
    close
      PRIMERINPUTFILE; ############ need to change these to $primerInputFH style
}

##############################
### Run blastn
if ($verbose) {
    print STDERR "starting BLAST\n";
}

# I did a bunch of testing and if you go over four processors for this step
# it actually slows things down.
# Also running ~100 fasta entries per batch seemed to be a good balance for RAM usage.
my $blastnProc;
if ( $processors <= 4 ) {
    $blastnProc = $processors;
}
else {
    $blastnProc = 4;
}

for ( my $i = 1 ; $i <= $primerOutFileNum ; $i++ ) {
    if ($verbose) {
        print STDERR "Running blast on file number " . $i . "\n";
    }

    my $blastCmd =
        $blastVer . " " . "-db " 
      . $blastDb . " "
      . "-query "
      . $tempDir . "/"
      . $tempName . "_"
      . $i . ".txt "
      . "-task blastn "
      . "-evalue 30000 "
      . "-word_size 7 "
      . "-max_hsps 100 "
      . "-max_target_seqs 50000 "
      . "-num_threads "
      . $blastnProc . " "
      . "-reward 1 "
      . "-penalty -1 "
      . "-gapopen 2 "
      . "-gapextend 1 "
      . "-outfmt \"6 qseqid sgi qlen qstart qend sstart send slen sstrand qseq sseq staxids\" "
      . "-ungapped "
      . "-sum_stats true";

    open BLASTRESULTS, "-|", $blastCmd or die "Blast query failed\n";

    while ( my $input = <BLASTRESULTS> ) {
        chomp $input;
        my @parsed = parseBlast($input);

        if ( scalar(@parsed) > 0 ) {
            my $shouldBeIncluded = 1;
            my (
                $qseqid,             $sgi,            $staxids,
                $ampStart,           $ampEnd,         $ampLength,
                $mismatchLoc,        $mismatch3Prime, $lastMismatchLoc,
                $lastMismatch3Prime, $qseq,           $sseq,
                $lastQseq,           $lastSseq
            ) = @parsed;
            if ( exists( $resultsHash{ $qseqid . "\t" . $sgi } ) ) {
                for my $result ( @{ $resultsHash{ $qseqid . "\t" . $sgi } } ) {
                    my @resultArray = split "\t", $result;

                    #if($resultArray[6] =~ /[fr]/) {
                    #    print STDERR join("\t", @resultArray), "\n";
                    #    die;
                    #}
                    if (   abs( $resultArray[3] - $ampStart ) < 20
                        && abs( $resultArray[4] - $ampEnd ) < 20 )
                    {
                        $shouldBeIncluded =
                          0;    # amplicon is the same as a previous one
                    }
                }
            }
            if ( $shouldBeIncluded == 1 ) {
                push @{ $resultsHash{ $qseqid . "\t" . $sgi } },
                  join( "\t", @parsed );
            }
        }
        if ( $verbose && $blastEntryCounter % 1000 == 0 ) {
            print STDERR
              "Parsing BLAST entry $blastEntryCounter                  \r";
        }
        $blastEntryCounter++;
    }
}

# print header
print join(
    "\t",
    "qSeqId",                    # query ID
    "sGi",                       # hit gi
    "sTaxids",                   # taxid of hit - can be multiple
    "sAmpStart",                 # start of amplicon location
    "sAmpEnd",                   # end amplicon location
    "sAmpLen",                   # amplicon length
    "primer1MismatchLoc",        # location of mismatches
    "primer1Mismatch3PrimeTip",  # number of mismatches in the 3' end
    "primer2MismatchLoc",        # location of mismatches other end
    "primer2Mismatch3Prime",     # number of mismatches in the 3' end, other end
    "primer1QSeq",               # query seq
    "primer1SSeq",               # hit seq
    "primer2QSeq",               # other end query seq
    "primer2SSeq"
  ),                             # other end hit seq
  "\n";

for my $primerGi ( keys %resultsHash ) {
    for my $result ( @{ $resultsHash{$primerGi} } ) {
        print $result, "\n";
    }
}

system( "rm " . $tempDir . "/" . $tempName . "_*.txt " ) if ( !$debug );

if ($verbose) {
    print STDERR "\nDone!\n";
}

print STDERR $numberOfBlastcmdCalls, " calls to blastdbcmd\n" if ($debug);

######################
### Subfunctions

sub addPrimersToHash {
    my $primer1 = $_[0];
    my $primer2 = $_[1];
    my $name    = $_[2];

    # Make up hash to use for parsing blast output
    if ( exists( $primerHash{$name} ) ) {
        die "Primer sets must have unique names\n";
    }
    else {
        @{ $primerHash{$name} } = ( $primer1, revComp($primer2) );
    }

    # deconvolute ambiguious bases for both primers
    my @primerFArray = ($primer1);
    my @primerRArray = ($primer2);

    if ( $primer1 =~ /[WSMKRYBDHVN]/ ) {
        @primerFArray = addDegeneratePrimer($primer1);
    }

    if ( $primer2 =~ /[WSMKRYBDHVN]/ ) {
        @primerRArray = addDegeneratePrimer($primer2);
    }

    my $primerCombinations  = 0;
    my $primerEntriesInFile = 0;
    my $primersOutFH;

    open $primersOutFH, ">",
      $tempDir . "/" . $tempName . "_" . $primerOutFileNum . ".txt";

    # Randomize the arrays in case I need to subsample
    # Don't shuffle if dubugging
    @primerFArray = shuffle(@primerFArray) if ( !$debug );
    @primerRArray = shuffle(@primerRArray) if ( !$debug );

    for my $seq1 (@primerFArray) {
        for my $seq2 (@primerRArray) {
            if ( $primerEntriesInFile >= $maxPrimersPerFile )
            {    # print out specified number of fasta entries per file
                $primerOutFileNum++;
                $primerEntriesInFile = 0;
                close $primersOutFH;
                open $primersOutFH, ">",
                  $tempName . "_" . $primerOutFileNum . ".txt";
            }

            if ( $primerCombinations < $maxPrimerVariantsTested )
            { # keep only the first $maxPrimerVariantsTested (2000) combinations
                print $primersOutFH ">", $name, "\n", $seq1, "N" x $nCount,
                  revComp($seq2), "\n";
                $primerEntriesInFile++;
            }
            $primerCombinations++;
        }
    }
    close $primersOutFH;

    if ($verbose) {
        print STDERR $primerCombinations, " primer combinations for ", $name,
          "\n";
        print STDERR "A maximum of ", $maxPrimerVariantsTested,
          " combinations will be tested.\n";
    }
}

sub addDegeneratePrimer {
    my $primer      = shift;
    my @primerArray = ($primer); #make an array containing the degenerate primer
    my @tempArray =
      ();    #make a temporary array to hold the new versions of the primers
    my $test = 1;
    while ( $test == 1 ) {
        for ( my $i = 0 ; $i < scalar(@primerArray) ; $i++ )
        {    # sort through primerArray
            if ( $primerArray[$i] =~ /[WSMKRYBDHVN]/ ) {
                push( @tempArray, getNewPrimerVersions( $primerArray[$i] ) );
            }
            else {
                push( @tempArray, $primerArray[$i] )
                  ;    #add normal primer to tempArray
            }
        }
        @primerArray = @tempArray;
        @tempArray   = ();
        if ( join( "", @primerArray ) =~ /[WSMKRYBDHVN]/ == 0 ) {
            $test = 0;
        }
    }
    return uniq(@primerArray);
}

sub getNewPrimerVersions {
    my $primer = shift;                    #primer sequence
    my @baseArray = split( "", $primer )
      ;    #split the primer up into individual bases in an array
    my @tempArray = ();
    for ( my $i = 0 ; $i < scalar(@baseArray) ; $i++ )
    {      #go through each base in the primer
        my $nucleotide = $baseArray[$i];
        if ( $nucleotide =~ /[WSMKRYBDHVN]/ )
        {    #if that base has a degenerate base
            for ( @{ $degeneratehash{$nucleotide} } )
            {    #go through the possible replacements
                my @copyArray =
                  @baseArray;    # copy this array so we can modify it
                $copyArray[$i] = $_
                  ; #then switch out the base in the copy array with a possibility
                push( @tempArray, join( "", @copyArray ) )
                  ; #and add the new decoded primer sequence to the end of the @primerArray
            }
        }
    }
    return @tempArray;
}

sub revComp {
    my $seq = shift;
    $seq =~ tr/ACGTacgtYRWSKMDVHBXN/TGCAtgcaRYWSMKHBDVXN/;
    $seq = reverse($seq);
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
    my (
        $qseqid, $sgi,  $qlen,    $qstart, $qend, $sstart,
        $send,   $slen, $sstrand, $qseq,   $sseq, $staxids
    ) = split "\t", $blastData;
    my @output;

    if ( $sgi eq $lastSgi ) {
        my (
            $lastQseqid,  $lastSgi,    $lastQlen, $lastQstart,
            $lastQend,    $lastSstart, $lastSend, $lastSlen,
            $lastSstrand, $lastQseq,   $lastSseq, $lastStaxids
        ) = split "\t", $lastLine;

        my @hitPositions =
          sort { $a <=> $b } ( $sstart, $send, $lastSstart, $lastSend );

        if (   $sstrand eq $lastSstrand
            && ( $hitPositions[3] - $hitPositions[0] <= $maxAmpLen )
            && ( $hitPositions[3] - $hitPositions[0] >= $minAmpLen ) )
        {

            # count primer mismatches
            my ( $lastMismatchLoc, $lastMismatch3Prime ) =
              mismatchCounter($lastLine)
              ; #$lastQseq, $lastSseq, $lastQstart, $lastQend, $lastQseqid, $lastSgi, $lastSstart);
            my ( $mismatchLoc, $mismatch3Prime ) = mismatchCounter($blastData)
              ;    #$qseq, $sseq, $qstart, $qend, $qseqid, $sgi);

            my $mismatchCountF = 0;
            if ( $lastMismatchLoc ne "" ) {
                $mismatchCountF = ( $lastMismatchLoc =~ tr/,/,/ );
                $mismatchCountF++;    # 3,4 only has one comma, so need to +1
            }
            my $mismatchCountR = 0;
            if ( $mismatchLoc ne "" ) {
                $mismatchCountR = ( $mismatchLoc =~ tr/,/,/ );
                $mismatchCountR++;
            }

            if (   $lastMismatch3Prime <= $primerTipMismatch
                && $mismatch3Prime <= $primerTipMismatch
                && $mismatchCountF <= $totalMismatchCount
                && $mismatchCountR <= $totalMismatchCount )
            {
                @output = (
                    $qseqid,             # query ID
                    $sgi,                # hit gi
                    $staxids,            # taxid of hit - can be multiple
                    $hitPositions[0],    # start of amplicon location
                    $hitPositions[3],    # end amplicon location
                    $hitPositions[3] - $hitPositions[0],    # amplicon length
                    $mismatchLoc,        # location of mismatches
                    $mismatch3Prime,     # number of mismatches in the 3' end
                    $lastMismatchLoc,    # location of mismatches other end
                    $lastMismatch3Prime
                    ,            # number of mismatches in the 3' end, other end
                    $qseq,       # query seq
                    $sseq,       # hit seq
                    $lastQseq,   # other end query seq
                    $lastSseq    # other end hit seq
                );
            }
        }
    }
    $lastLine = $blastData;
    $lastSgi  = $sgi;
    return @output;
}

sub mismatchCounter {
    my $blastData = shift;
    my (
        $qseqid, $sgi,  $qlen,    $qstart, $qend, $sstart,
        $send,   $slen, $sstrand, $qseq,   $sseq, $staxids
    ) = split "\t", $blastData;

    my $primerShouldStart;
    my $primerShouldEnd;
    my $primerArrayIndex = 0;
    my $primerDir        = "for:";
    my $primerRealSeq;
    my $hitRealSeq     = $sseq;
    my $mismatchLocs   = "";
    my $mismatch3Prime = 0;

# Check to be sure both sequence and primer are full length and fill in if needed

    if ( $qstart < 20 ) {    # query is the forward primer
        $primerShouldStart = 1;
        $primerShouldEnd   = length( @{ $primerHash{$qseqid} }[0] );
        $primerRealSeq     = $forward;
        if ( $sstrand eq "minus" ) {
            $primerRealSeq = revComp($primerRealSeq);
        }
    }
    else {                   # query is the reverse primer
        $primerShouldStart =
          length( @{ $primerHash{$qseqid} }[0] ) + $nCount + 1;
        $primerShouldEnd =
          $primerShouldStart + length( @{ $primerHash{$qseqid} }[1] ) - 1;
        $primerArrayIndex = 1;
        $primerDir        = "rev:";
        $primerRealSeq    = $reverse;

        if ( $sstrand eq "plus" ) {
            $primerRealSeq = revComp($primerRealSeq);
        }

    }

    # Check left end and add sequences as needed
    print STDERR $blastData, "\n", $primerHash{$qseqid}[$primerArrayIndex],
      " real primer\n", $primerRealSeq, "\t", $hitRealSeq,
      "primerBeforeCheckingStart\n"
      if ($debug);
    if ( $qstart != $primerShouldStart ) {
        if ( $qstart == $primerShouldStart + 1 ) {
            $hitRealSeq = "Z" . $hitRealSeq;
        }
        else {
            my $seqBaseNumberToGet = $sstart -
              ( $qstart - $primerShouldStart );    #base position is 1-based
            if ( $seqBaseNumberToGet > 0 ) {
                my $seqRange = $seqBaseNumberToGet . "-" . ( $sstart - 1 );

                my $getSeqCmd =
                    "blastdbcmd -target_only -outfmt \"%s\" -db " 
                  . $blastDb
                  . " -entry "
                  . $sgi
                  . " -range "
                  . $seqRange;
                print STDERR $getSeqCmd, "For start\n" if ($debug);
                my $seqBases = `$getSeqCmd`;
                $numberOfBlastcmdCalls++;
                chomp $seqBases;
                $seqBases   = uc($seqBases);
                $hitRealSeq = $seqBases . $hitRealSeq;
            }
            elsif ( ( $sstart - 1 ) > 0 ) {

                # get what I can and put on N's for the rest
                my $seqRange = "1-" . ( $sstart - 1 );

                my $getSeqCmd =
                    "blastdbcmd -target_only -outfmt \"%s\" -db " 
                  . $blastDb
                  . " -entry "
                  . $sgi
                  . " -range "
                  . $seqRange;
                print STDERR $getSeqCmd, "For partial start\n" if ($debug);
                my $seqBases = `$getSeqCmd`;
                $numberOfBlastcmdCalls++;
                chomp $seqBases;
                $seqBases   = uc($seqBases);
                $hitRealSeq = $seqBases . $hitRealSeq;

                # Add N's for the rest
                $hitRealSeq =
                  "N" x ( abs($seqBaseNumberToGet) + 1 ) . $hitRealSeq;
            }
            else {

       # the missing sequence from the subject is just missing from the sequence
                $hitRealSeq =
                  "N" x ( abs($seqBaseNumberToGet) + 1 ) . $hitRealSeq;
            }
        }
    }

    # Check right end of primer and add sequences as needed
    print STDERR $primerRealSeq, "\t", $hitRealSeq, "primer after start fix\n"
      if ($debug);
    if ( $qend != $primerShouldEnd ) {
        if ( $qend == $primerShouldEnd - 1 ) {
            $hitRealSeq = $hitRealSeq . "Z";
        }
        else {
            my $seqEndNumberToGet =
              $send + ( $primerShouldEnd - $qend );    #base position is 1-based
            my $seqRange = ( $send + 1 ) . "-" . $seqEndNumberToGet;

            my $getSeqCmd =
                "blastdbcmd -target_only -outfmt \"%s\" -db " 
              . $blastDb
              . " -entry "
              . $sgi
              . " -range "
              . $seqRange;
            my $seqBases = `$getSeqCmd`;
            $numberOfBlastcmdCalls++;
            print STDERR $getSeqCmd, "\tcommand to get end\n" if ($debug);
            chomp $seqBases;
            $seqBases   = uc($seqBases);
            $hitRealSeq = $hitRealSeq . $seqBases;

  # Check if the full sequence was not available, and fill in with N's otherwise
            if ( length($hitRealSeq) < length($primerRealSeq) ) {
                $hitRealSeq =
                  "N" x ( length($primerRealSeq) - length($hitRealSeq) )
                  . $hitRealSeq;
            }
        }
    }

    print STDERR $primerRealSeq, "\t", $hitRealSeq, "\t", $primerDir,
      " after end fix\n"
      if ($debug);

    # Count the mismatches in the aligned portion of the primer
    # no need to count mismatches if identical
    # This doesn't work if ambiguous bases are present, but that's fine
    if ( $primerRealSeq ne $hitRealSeq ) {

        # if the hit is on the negative strand, need to revcomp the sequence
        if ( $sstrand ne "plus" ) {
            my $primerRealSeq = revComp($primerRealSeq);
            my $hitRealSeq    = revComp($hitRealSeq);
        }

# rev complement forward primers so they're facing 3'---5' so $i = 0 is the 3' end
        if ( $primerDir eq "for:" ) {
            $primerRealSeq = revComp($primerRealSeq);
            $hitRealSeq    = revComp($hitRealSeq);
        }

        # Do the comparison for each base
        for ( my $i = 0 ; $i < length($primerRealSeq) ; $i++ ) {
            my $qBase = substr( $primerRealSeq, $i, 1 );
            my $sBase = substr( $hitRealSeq,    $i, 1 );
            if ( $sBase !~ /[$degenerateRegexHash{$qBase}]/ ) {

         # Need to complement the mismatch since we revComp'd the sequence above
                $mismatchLocs .= "," . ( $i + 1 ) . "_" . revComp($sBase);
                if ( $i <= $primerTipLen ) {
                    $mismatch3Prime++;
                }
            }
        }
    }
    $mismatchLocs =~ s/^,//;    # get rid of leading comma
    return ( $primerDir . $mismatchLocs, $mismatch3Prime );
}

##############################
# POD
##############################

#=pod

=head SYNOPSIS

Summary:

    xxxxxx.pl - description

Usage:

    perl xxxxxx.pl [options]


=head OPTIONS

Options:

    --verbose
    --help

=cut
