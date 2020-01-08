#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use DBI;


##############################
# By Matt Cannon
# Date: a while ago
# Last modified: 8-19-19
# Title: bsPrimerTree.pl
# Purpose: local version (sortof) of the PrimerTree R package
##############################


##############################
### To do maybe
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
my $taxDb;
my $maxAlignedSeqs = "inf";
my $maxSeqPerSp = "inf";
my $noPlots;
my $plotPtSz = 12;
my $plotFtSz = 12;
my $minTaxaToPlot = 20;

# i = integer, s = string
GetOptions ("verbose"             => \$verbose,
            "help"                => \$help,
	    "threads=i"           => \$threads,
	    "inFile=s"            => \$inFile,
	    "outDir=s"            => \$outDir,
	    "blastDb=s"           => \$blastDb,
	    "taxDb=s"             => \$taxDb,
	    "maxAlignedSeqs=i"    => \$maxAlignedSeqs,
	    "maxSeqsPerSpecies=i" => \$maxSeqPerSp,
	    "noPlots"             => \$noPlots,
	    "plotPointSize=i"     => \$plotPtSz,
	    "plotFontSize=i"      => \$plotFtSz,
	    "maxTaxaToPlot=i"     => \$minTaxaToPlot            
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %seqTrimHash;
my %taxNamesHash;
my %taxCountHash;
my %mismatchHash;
my %speciesSeqCountHash;

##############################
# Code
##############################

if(-w "." eq "") { # check for write permissions
      print STDERR "No write permissions to this directory!\n\n";
      die;
}
######### Put other input checks here ##########
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
    my @time = localtime(time);
    print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
}

open my $blastDbInputFile, ">", $outDir . "seqsToGet.txt";
open my $ampliconLenFH, ">", $outDir . "ampliconLengths.txt";

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
       $primer2Mismatch5PrimeTip,              # number of mismatches in the 5' end, other end
       $primer1QSeq,                           # query seq
       $primer1SSeq,                           # hit seq
       $primer2QSeq,                           # other end query seq
       $primer2SSeq)                           # other end hit seq
	= split "\t", $input;

    print $blastDbInputFile $sGi, "\t", $sAmpStart, "-", $sAmpEnd, "\n";
    print $ampliconLenFH $sGi, "\t", $sAmpLen, "\n";
    
    # Store mismatch info to parse out later once we have the taxa info
    $primer1MismatchCount =~ s/:/\t/;
    $primer2MismatchCount =~ s/:/\t/;
    

    $mismatchHash{$sTaxids}{$primer1MismatchCount . "\t" . $primer1Mismatch5PrimeTip}++; 
    $mismatchHash{$sTaxids}{$primer2MismatchCount . "\t" . $primer2Mismatch5PrimeTip}++;

    # Store hit seq and other end hit seqs for each hit so I can trim them off later
    @{ $seqTrimHash{$sGi} } = ($primer1SSeq, $primer2SSeq);
}

close $inputFH;
close $blastDbInputFile;

######################
### Get sequence for each hit
### Also get taxonomic info and put into header

if($verbose) {
    print STDERR "Getting amplicon sequence information using blastdbcmd.\n";
    my @time = localtime(time);
    print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
}

### connect to taxonomyDb created using makeTaxonomyDb.pl
my $dsn      = "dbi:SQLite:dbname=$taxDb";
my $user     = "";
my $password = "";
my $dbh = DBI->connect($dsn, $user, $password, {
   PrintError       => 0,
   RaiseError       => 1,
   AutoCommit       => 1,
});

### Set up taxa db query
my $query = 'SELECT species, genus, family, "order", 
		class, phylum, kingdom, superkingdom, tax_name FROM 
		taxonomy WHERE tax_id == ?';
my $sth = $dbh->prepare($query);

### Put together command to get sequences 
### output will have header with >gi_taxid \n sequence for each query
### perl code switches "@" with newline
my $blastDBCmdCmd = "blastdbcmd " . 
    "-target_only " . 
    "-db " . $blastDb . 
    " -entry_batch " . $outDir . "seqsToGet.txt " . 
    "-outfmt \">\%g_\%T\@\%s\" " .
    "| perl -pe \'s/\@/\n/\' "; 

### Array of data to print to speed up the program
my @printable;

### Go through sequences retrieved from blastDb and get the taxa info as it goes
$/ = "\n>";
open my $fastaWithTaxaFile, ">", $outDir . "seqsWithTaxa.fasta";
open my $blastDbCmdResponse, "-|", $blastDBCmdCmd or die "blastdbcmd query failed\n";

open my $mismatchFile, ">", $outDir . "primerMismatches.txt";
print $mismatchFile join("\t", "taxid", "superkingdom", "kingdom", "phylum", 
			 "class", "order", "family", "genus", "species", 
			 "direction", "mismatchTotal", "mismatch5Prime"), "\n";


while(my $blastInput = <$blastDbCmdResponse>) {
    chomp $blastInput;
    $blastInput =~ s/^>//;
    my ($header, $seq) = split "\n", $blastInput;

    my $gi = $header;
    $gi =~ s/_.+//;
    
    my $taxid = $header;
    $taxid =~ s/.+_//;

    ### use the subject sequences given by bsPrimerBlast to trim primer sequence aligned portion off sequence 

    ###
    ############  See Haikel's Nege2 primers, but the returned sequence/primers may need to be rev comp'd to find the match for trimming
    ### Negevirus_group2_F GTTGCWGGTCACGGTAARAC	Negevirus_group2_R CRTCAGCWGGAATWCGATAC

    ### Check this if else stack for logic later ############################

    my @seqsToTrimOff = @{ $seqTrimHash{$gi} };

    my @revCompSeqsToTrimOff = (revComp($seqsToTrimOff[0]), revComp($seqsToTrimOff[1])); 

    ###### I may need to just align the primers to the sequence instead of including all possible orientations here.... 
    if($seq =~ /^$seqsToTrimOff[0].+$seqsToTrimOff[1]$/) { 
	$seq =~ s/^$seqsToTrimOff[0]//;
	$seq =~ s/$seqsToTrimOff[1]$//;
    } elsif($seq =~ /^$revCompSeqsToTrimOff[1].+$revCompSeqsToTrimOff[0]$/) {
	$seq =~ s/^$revCompSeqsToTrimOff[1]//;
	$seq =~ s/$revCompSeqsToTrimOff[0]$//;
    } elsif($seq =~ /^$seqsToTrimOff[1].+$seqsToTrimOff[0]$/) {
	$seq =~ s/^$seqsToTrimOff[1]//;
	$seq =~ s/$seqsToTrimOff[0]$//;
    } elsif($seq =~ /^$revCompSeqsToTrimOff[0].+$revCompSeqsToTrimOff[1]$/) {
	$seq =~ s/^$revCompSeqsToTrimOff[0]//;
	$seq =~ s/$revCompSeqsToTrimOff[1]$//;
    } else {
	print STDERR "Warning! Subject sequences don't match sequence returned by blastdbcmd $header\t", join(",", @seqsToTrimOff), "\t", $seq, "\n"
    }

    ### Get taxa info from database
    $sth->execute($taxid); 
    
    my ($species, $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $tax_name);
    while(my $row = $sth->fetchrow_hashref){
	$species      = "$row->{species}";
	$genus        = "$row->{genus}";
	$family       = "$row->{family}";
	$order        = "$row->{order}";
	$class        = "$row->{class}";
	$phylum       = "$row->{phylum}";
	$kingdom      = "$row->{kingdom}";
	$superkingdom = "$row->{superkingdom}";
	$tax_name     = "$row->{tax_name}";
    }


    ### print out 
    if(defined($species)) {
	# some of the entries in the tax db don't have species labeled >:-|
	if($species eq "NA") {
	    $species = $tax_name;
	}
	$speciesSeqCountHash{$species}++;
	
	my $newHeader = ">"   . # header with taxa info and unique number to make alignment program happy
			"s-"  . $species . ":" .
			"sk-" . $superkingdom . ":" .
			"k-"  . $kingdom . ":" .
			"p-"  . $phylum . ":" .
			"c-"  . $class . ":" .
			"o-"  . $order . ":" .
			"f-"  . $family . ":" .
			"g-"  . $genus . ":" . 
			"_" . $speciesSeqCountHash{$species};
			
	# put it into printable array if I haven't seen too many of that species
	if($speciesSeqCountHash{$species} < $maxSeqPerSp) {
	    push @printable, $newHeader . "\n" . $seq;
	}
	
	# record how many of each taxa I've seen to print out later
	$taxCountHash{$superkingdom . "\t" . 
			  $kingdom . "\t" . 
			  $phylum . "\t" . 
			  $class . "\t" . 
			  $order . "\t" . 
			  $family . "\t" . 
			  $genus . "\t" . 
			  $species}++;
	
	$taxNamesHash{species}{$species}++;
	$taxNamesHash{genus}{$genus}++;
	$taxNamesHash{family}{$family}++;
	$taxNamesHash{order}{$order}++;
	$taxNamesHash{class}{$class}++;
	$taxNamesHash{phylum}{$phylum}++;
	$taxNamesHash{kingdom}{$kingdom}++;
	$taxNamesHash{superkingdom}{$superkingdom}++;


	# Print out information on how many mismatches there are
	for my $mismatch (keys %{ $mismatchHash{$taxid} } ) {
	    print $mismatchFile $taxid . "\t" . 
		                $superkingdom . "\t" . 
				$kingdom . "\t" . 
				$phylum . "\t" . 
				$class . "\t" . 
				$order . "\t" . 
				$family . "\t" . 
				$genus . "\t" . 
				$species . "\t" .
				$mismatch . "\n"; 
	}

    } else {
	print STDERR "Taxonomy not found for gi: $gi, taxid: $taxid\n";
	print STDERR "Make sure your blast database and taxonomy database " . 
			"were downloaded at about the same time\n";
    }
}

close $mismatchFile;
close $blastDbCmdResponse;
$dbh->disconnect;


### Don't want to print out too many sequences for the alignment, so print out randomly if there are too many
if($maxAlignedSeqs < scalar(@printable)) {
    # random print
    my $printChance = $maxAlignedSeqs / scalar(@printable);
    for my $entry (@printable) {
	if(rand() <= $printChance) {
	    print $fastaWithTaxaFile $entry, "\n";
	} 
    }
    # otherwise just print everything
} else {
    print $fastaWithTaxaFile join("\n", @printable), "\n";
}


close $fastaWithTaxaFile;

        ##### Need to keep in mind that the output fasta may  have multiple ">"s in the header
$/ = "\n";



########################
### Print out a log of the number of taxa within each taxonomic level
open my $taxaSummaryFile, ">", $outDir . "taxaSummary.txt";
print $taxaSummaryFile "superkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tCount\n";

#for my $level ("subkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species") {
#    for my $taxName (keys %{ $taxNamesHash{$level} } ) {
#	print $taxaSummaryFile $level, "\t", $taxName, "\t", $taxNamesHash{$level}{$taxName}, "\n";
for my $taxName (keys %taxCountHash) {
    print $taxaSummaryFile $taxName, "\t", $taxCountHash{$taxName}, "\n";
    
}
#}

close $taxaSummaryFile;

########################
### Align reads
if($verbose) {
    print STDERR "Aligning reads with Mafft.\n";
    my @time = localtime(time);
    print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
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
### Make tree with FastTree
if(!$noPlots) {
    if($verbose) {
	print STDERR "Making tree with FastTree.\n";
	my @time = localtime(time);
	print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
    }
    
    my $treeCmd = "FastTree -quote -quiet -nt " . $outDir . "seqsWithTaxaAligned.fasta > " . $outDir . "tree.nwk";
    
    if($verbose) {
	$treeCmd =~ s/-quiet //;
    }
    
    system($treeCmd);
    
    
    #######################
    ### Make Dendroscope instructions
    if($verbose) {
	print STDERR "Making Dendroscope command files and running Dendroscope\n";
	print STDERR "Dendroscope may open and do random stuff. Do not close it or interact with the window or the fabric of time will come unraveled\n";
	my @time = localtime(time);
	print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
    }
    
    for my $taxaLevel (keys %taxNamesHash) {
	my %taxaLevelNameHash = (
	    "species" 	=> "s",
	    "superkingdom"	=> "sk",
	    "kingdom" 	=> "k",
	    "phylum"	=> "p",
	    "class"		=> "c",
	    "order"		=> "o",
	    "family"	=> "f",
	    "genus"		=> "g");
	
	my $dendroHeader = join("\n", 
				"open file=" . $outDir . "tree.nwk" . ";",
				"set window width=3000 height=3000;", # Have to set this to a high value to avoid collapsed branches :-(
				"go tree=1;", 
				"show nodelabels=false;",
				"show edgelabels=false;",
				"set drawer=RadialPhylogram;",
				"set scale=1000;",
				"deselect all;",
				#"replace searchtext=^R replacetext=_ all=true regex=true;",
				"replace searchtext=_\\d+\$ replacetext=_ all=true regex=true;",
				"select nodes=leaves;",
				"set nodeshape=oval;" ,
				"set nodesize=" . $plotPtSz . ";" ,
				"set fontsize=" . $plotFtSz . ";") . "\n";
	
	my $dendroTail = join("\n", 
			      "set radiallabels=true;",
			      "show scalebar=true;",
			      "deselect all;",
			      "exportimage file=\'" . $outDir . "treePlot_" . $taxaLevel . ".svg\' textasshapes=false format=SVG replace=true;",
			      "quit;") . "\n";
	
	open my $dendroInstFile, ">", $outDir . "dendroInstructionFile_" . $taxaLevel . ".txt";
	
	my @taxaList = keys %{ $taxNamesHash{$taxaLevel} }; 
	
	my $numberOfNames = scalar(@taxaList);
	my $colorSteps = int($numberOfNames ** (1/3)) + 1;
	my $colMaxVal = 255;
	my $stepSize = int($colMaxVal / $colorSteps);
	my @colors;
	
	#modify dendroHeader to cut everything except the target taxa level. - Move this to $dendroHeader declaration?
	$dendroHeader .= 
	    "replace searchtext=[Rs].*" . $taxaLevelNameHash{$taxaLevel} . "- replacetext=---------- all=true regex=true;\n" .
	    "replace searchtext=:.* replacetext=------------ all=true regex=true;\n";
	
	# loop through to make all combos
	for(my $i = 255; $i > 0; $i -= $stepSize){
	    for(my $j= 255; $j > 0; $j -= $stepSize    ) {
		for(my $k= 255; $k > 0; $k -= $stepSize     ){
		    if($i != 255 || $j != 255 || $k != 255) { # don't want white
			push @colors, $i . " " . $j . " " . $k;
		    }
		}
	    }
	}
	
	print $dendroInstFile $dendroHeader;
	
	for my $taxName (@taxaList) {
	    if($taxNamesHash{$taxaLevel}{$taxName} > $minTaxaToPlot) {
		# should I keep only 10% or so of the labels?
		my $colToUse = shift @colors;
		
		print $dendroInstFile "find searchtext=\'" . $taxName, "\';\n";
		print $dendroInstFile "show nodelabels=true;";
		print $dendroInstFile "set color=" . $colToUse . ";\n";
		print $dendroInstFile "set fillcolor=" . $colToUse . ";\n";
		print $dendroInstFile "set labelcolor=" . $colToUse . ";\n";
		print $dendroInstFile "deselect all;\n";
	    } 
	}
	
	print $dendroInstFile $dendroTail;
	close $dendroInstFile;
	my $dendroscopeCmd = "Dendroscope -g -c " . $outDir . "dendroInstructionFile_" . $taxaLevel . ".txt";
	if($verbose) {
	    print STDERR "starting Dendroscope for $taxaLevel\n";
	}
	my $response = `$dendroscopeCmd`; 
	system($dendroscopeCmd);
	if($verbose) {
	    print STDERR "Done with Dendroscope\n";
	}
	# Uncomment when a working version of ImageMagick is available. Or don't. Whatever.
	#my $imgConvertCmd = "convert -density 300 " . $outDir . "treePlot_" . $taxaLevel . ".svg " . $outDir . "treePlot_" . $taxaLevel . ".jpg";
	#system($imgConvertCmd);
}



# check if imageMagick is installed -> run convert --version and grep for ImageMagick
# convert -density 300 testFigure.svg testConvert.jpg

}

if($verbose) {
    print STDERR "Done!!!!\n";
    my @time = localtime(time);
    print STDERR "Time: ", $time[2] . ":" . $time[1], "\n";
}
    

sub revComp {
    my $seq = shift;
    $seq =~ tr/ACGTacgtYRWSKMDVHBXN/TGCAtgcaRYWSMKHBDVXN/;
    $seq = reverse ($seq);
    return $seq;
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
