#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

use DBI;


##############################
# By Matt Cannon
# Date: 
# Last modified: 
# Title: .pl
# Purpose: 
##############################


##############################
### To do maybe
# limit number of hits?

# limit taxonomically using new version of blast?

# implement checks on input and successful run of each command. 


##############################


##############################
# Options
##############################


my $verbose;
my $help;
my $threads = 1;
my $inFile;
my $outDir;
my $blastDb;
my $maxSeqPerSp = "inf";
my $plotPtSz = 4;
my $plotFtSz = 6;

# i = integer, s = string
GetOptions ("verbose"             => \$verbose,
            "help"                => \$help,
	    "threads=i"           => \$threads,
	    "inFile=s"            => \$inFile,
	    "outDir=s"            => \$outDir,
	    "blastDb=s"           => \$blastDb,
	    "maxSeqsPerSpecies=i" => \$maxSeqPerSp,
	    "plotPointSize=i"     => \$plotPtSz,
	    "plotFontSize=i"      => \$plotFtSz
            
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %taxidHash;
my %taxNamesHash;
my %fastaHash;

##############################
# Code
##############################

if(-w "." eq "") { # check for write permissions
      print STDERR "No write permissions to this directory!\n\n";
      die;
}

if($outDir eq "") {
    print STDERR "please provide an output directory\n\n";
    die;
} else {
      $outDir .= "/";
}

mkdir $outDir;

if($verbose) {
    print STDERR "Output directory: ", $outDir, "\n\n";
}


##############################
### Read in bsPrimerBlast file and write out files for getTaxa and to
### to get the sequences using blastdbcmd
if($verbose) {
    print STDERR "Parsing bsPrimerBlast.pl input and preparing files for analysis.\n";
}

open my $blastDbInputFile, ">", $outDir . "seqsToGet.txt";

open my $inputFH, "$inFile" or die "Could not open input\nWell, crap\n";
my $header = <$inputFH>; 

while (my $input = <$inputFH>){
    chomp $input;
    my($qSeqId,                                # query ID
       $sGi,                                   # hit gi
       $sTaxids,                               # taxid of hit - can be multiple
       $sSciNames,                             # scientific name
       $sComNames,                             # common name
       $sAmpStart,                             # start of amplicon location
       $sAmpEnd,                               # end amplicon location
       $sAmpLen,                               # amplicon length
       $primer1MismatchCount,                  # number of mismatches
       $primer1Mismatch5PrimeTip,              # number of mismatches in the 5' end
       $primer2MismatchCount,                  # number of mismatches other end
       $primer2Mismatch5Prime,                 # number of mismatches in the 5' end, other end
       $primer1QSeq,                           # query seq
       $primer1SSeq,                           # hit seq
       $primer2QSeq,                           # other end query seq
       $primer2SSeq)                           # other end hit seq
	= split "\t", $input;

    print $blastDbInputFile $sGi, "\t", $sAmpStart, "-", $sAmpEnd, "\n";
    
    my @taxa = split ";", $sTaxids;
    for(@taxa) {
	$taxidHash{$_} = 1;
    }
}

close $inputFH;
close $blastDbInputFile;



######################
###  get sequence for each hit
if($verbose) {
    print STDERR "Getting amplicon sequence information using blastdbcmd.\n";
}

my $blastDBCmdCmd = "blastdbcmd " . 
    "-db " . $blastDb . 
    " -entry_batch " . $outDir . "seqsToGet.txt " . 
    "-outfmt \">\%T\@\%s\" " .
    "| perl -pe \'s/\@/\n/\' "; # .
#                        "> " . $outDir . "seqs.fasta";
    
#                        system($blastDBCmdCmd); 
$/ = "\n>";
open my $blastDbCmdResponse, "-|", $blastDBCmdCmd or die "blastdbcmd query failed\n";
while(my $blastInput = <$blastDbCmdResponse>) {
   chomp $blastInput;
    $blastInput =~ s/^>//;
    my ($header, $seq) = split "\n", $blastInput;
    $fastaHash{$header}{$seq} = 1;
}


        ##### Need to keep in mind that the output fasta may  have multiple ">"s in the header
$/ = "\n";

######################
### get taxa info for each hit
if($verbose) {
    print STDERR "Getting taxonomic information using getTaxa.pl.\n";
}


my @taxaToGet = keys %taxidHash;

open my $getTaxaInputFile, ">", $outDir . "taxaToGet.txt";
print $getTaxaInputFile join("\n", @taxaToGet), "\n";


my $getTaxaCmd = "perl ~/bin/getTaxa.pl --taxid " . $outDir . "taxaToGet.txt"; 

open my $getTaxaResponse, "-|", $getTaxaCmd or die "getTaxa.pl query failed\n";

my $taxaHeader = <$getTaxaResponse>;
my @taxaHeaderArray = split "\t", $taxaHeader;
@taxaHeaderArray = map(substr($_, 0, 1), @taxaHeaderArray);


#######################
### Combine fasta sequences and taxa info to write out new fasta
if($verbose) {
    print STDERR "Combining taxonomic information with fasta sequences.\n";
}


open my $fastaWithTaxaFile, ">", $outDir . "seqsWithTaxa.fasta";

while(my $input = <$getTaxaResponse>) {
    chomp $input;
    my @taxonomyArray = split "\t", $input;

    my $newHeader = ">";
    for(my $i = 1; $i < scalar(@taxonomyArray); $i++) {
	$newHeader .= $taxaHeaderArray[$i] . "-" . $taxonomyArray[$i] . ":";
	$taxNamesHash{$taxaHeaderArray[$i]}{$taxonomyArray[$i]}++;       # Make a hash of all the taxonomy names found at each level
    }
    
    my $seqCount = 1;
    for my $seq (keys %{ $fastaHash{$taxonomyArray[0]} }) {
	if($seqCount < $maxSeqPerSp) {
	    print $fastaWithTaxaFile $newHeader, "\n", $seq, "\n";
	}
	$seqCount++;
    }
}

# Need to make instructions for dendroscope here using taxa info


######### I should print out a log of the number of taxa within each taxonomic level, or maybe just a table of all taxa?


########################
### Align reads
if($verbose) {
    print STDERR "Aligning reads with Mafft.\n";
}

my $mafftCmd = "mafft --quiet --auto --adjustdirectionaccurately " . 
    "--thread " . $threads . " " .
    $outDir . "seqsWithTaxa.fasta > " . 
    $outDir . "seqsWithTaxaAligned.fasta";

if($verbose) {
    $mafftCmd =~ s/--quiet //;
}

system($mafftCmd);


#######################
### Make tree with mega_cc

if($verbose) {
    print STDERR "Making tree with Mega.\n";
}

# Add number of threads to mao file
my $maoSwitchCmd = "perl -pe \'s/@/" . $threads . "/\' megaCmdsStub.mao > megaCmds.mao";
system($maoSwitchCmd);

my $megaCmd = "megacc -s -a megaCmds.mao -d " . $outDir . "seqsWithTaxaAligned.fasta -o " . $outDir . "tree.nwk";

if($verbose) {
    $megaCmd =~ s/-s //;
}

###### Check if the tree.nwk file already exists, and if it does, delete it to avoid mega screwing up via tree(2).nwk

system($megaCmd);


#######################
### Make Dendroscope instructions
if($verbose) {
    print STDERR "Making Dendroscope command files and running Dendroscope\n";
    print STDERR "Dendroscope may open and do random stuff. Do not close it or interact with the window or the fabric of time will come unraveled\n";
}

for my $taxaLevel (keys %taxNamesHash) {
    my $dendroHeader = join("\n", 
			    "open file=" . $outDir . "tree.nwk" . ";",
			    "set window width=2000 height=2000;",
			    "go tree=1;", 
			    "set drawer=RadialPhylogram\;",
			    "replace searchtext=^R replacetext=_ all=true regex=true;",
			    "select nodes=leaves;",
			    "set nodeshape=oval;" ,
			    "set nodesize=" . $plotPtSz . ";" ,
			    "set fontsize=" . $plotFtSz . ";") . "\n";
    
    my $dendroTail = join("\n", 
			  "deselect all;",
			  "set radiallabels=true;",
			  "show scalebar=true;",
			  "exportimage file=\'" . $outDir . "treePlot_" . $taxaLevel . ".svg\' textasshapes=false format=SVG replace=true;",
			  "quit;") . "\n";

    

    open my $dendroInstFile, ">", $outDir . "dendroInstructionFile_" . $taxaLevel . ".txt";
    
    my @taxaList = keys %{ $taxNamesHash{$taxaLevel} }; 

    my $numberOfNames = scalar(@taxaList);
    my $colorSteps = int($numberOfNames ** (1/3)) + 1;
    my $colMaxVal = 255;
    my $stepSize = int($colMaxVal / $colorSteps);
    my @colors;

    #modify dendroHeader to cut everything except the target taxa level
    $dendroHeader .= 
	"replace searchtext=[Rs].*" . $taxaLevel . "- replacetext=- all=true regex=true;\n" .
	"replace searchtext=:.* replacetext=- all=true regex=true;\n";

    # loop through to make all combos
    for(my $i = 0; $i < 255; $i += $stepSize){
	for(my $j= 0; $j < 255; $j += $stepSize    ) {
	    for(my $k= 0; $k < 255; $k += $stepSize     ){
		if($i != 0 || $j != 0 || $k != 0) { # don't want pure white
		    push @colors, $i . " " . $j . " " . $k;
		}
	    }
	}
    }

    print $dendroInstFile $dendroHeader;
    for my $taxName (@taxaList) {
	# should I keep only 10% or so of the labels?
	my $colToUse = shift @colors;

	print $dendroInstFile "find searchtext=\'" . $taxName, "\';\n";
	print $dendroInstFile "set color=" . $colToUse . ";\n";
	print $dendroInstFile "set fillcolor=" . $colToUse . ";\n";
	print $dendroInstFile "set labelcolor=" . $colToUse . ";\n";
	print $dendroInstFile "deselect all;\n";
    } 


    print $dendroInstFile $dendroTail;
    close $dendroInstFile;
    my $dendroscopeCmd = "Dendroscope -g -c " . $outDir . "dendroInstructionFile_" . $taxaLevel . ".txt";
    system($dendroscopeCmd);
    my $imgConvertCmd = "convert -density 300 " . $outDir . "treePlot_" . $taxaLevel . ".svg " . $outDir . "treePlot_" . $taxaLevel . ".jpg";
    system($imgConvertCmd);
}



# check if imageMagick is installed -> run convert --version and grep for ImageMagick
# convert -density 300 testFigure.svg testConvert.jpg

if($verbose) {
    print STDERR "Done!!!!\n";
}

##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    bsPrimerTree.pl - knockoff of the R package primerTree
    
Usage:

    perl bsPrimerTree.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
