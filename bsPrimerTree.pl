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
my %mismatchLocsHash;
my %speciesSeqCountHash;
my %alignedSeqHash;
my %distHash;
my %ampLenHash;
my %taxaHash;
my %giTaxaHash;
my $blastHitCount = 0;


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
       $primer1MismatchLocs,                   # locations of mismatches
       $primer1Mismatch3PrimeTip,              # number of mismatches in the 3' end
       $primer2MismatchLocs,                   # locations of mismatches other end
       $primer2Mismatch3PrimeTip,              # number of mismatches in the 3' end, other end
       $primer1QSeq,                           # query seq
       $primer1SSeq,                           # hit seq
       $primer2QSeq,                           # other end query seq
       $primer2SSeq)                           # other end hit seq
	= split "\t", $input;

    print $blastDbInputFile $sGi, "\t", $sAmpStart, "-", $sAmpEnd, "\n";

    # some taxids have two taxids separated by ";" - need to deal with this better at some point 
    $sTaxids =~ s/;.+//; 
    
    # store amplicon length until I have taxa info
    $ampLenHash{$sTaxids}{$sAmpLen}++;
    
    # Store mismatch info to parse out later once we have the taxa info
    # Primer 1 mismatches
    my ($primer1Dir, $primer1MismatchLocNums) = split ":", $primer1MismatchLocs;
    my @primer1MismatchLocArray = split ",", $primer1MismatchLocNums;
    my $primer1MismatchCount = scalar(@primer1MismatchLocArray);

    $mismatchHash{$sTaxids}{$primer1Dir . "\t" . $primer1MismatchCount . "\t" . $primer1Mismatch3PrimeTip}++;
    $mismatchLocsHash{$sTaxids}{$primer1Dir . "\t" . join(",", @primer1MismatchLocArray)}++;

    # Primer 2 mismatches
    my ($primer2Dir, $primer2MismatchLocNums) = split ":", $primer2MismatchLocs;
    my @primer2MismatchLocArray = split ",", $primer2MismatchLocNums;
    my $primer2MismatchCount = scalar(@primer2MismatchLocArray);
    
    $mismatchHash{$sTaxids}{$primer2Dir . "\t" . $primer2MismatchCount . "\t" . $primer2Mismatch3PrimeTip}++;
    $mismatchLocsHash{$sTaxids}{$primer2Dir . "\t" . join(",", @primer2MismatchLocArray)}++;


    # Store hit seq and other end hit seqs for each hit so I can trim them off later
    @{ $seqTrimHash{$sGi} } = ($primer1SSeq, $primer2SSeq);
    
    #$blastHitCount++; # store this for primer mismatch printing
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
my $query = 'SELECT species_level, genus_level, family_level, order_level, 
		class_level, phylum_level, kingdom_level, superkingdom_level, tax_name FROM 
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

while(my $blastInput = <$blastDbCmdResponse>) {
    chomp $blastInput;
    $blastInput =~ s/^>//;
    my ($header, $seq) = split "\n", $blastInput;

    my $gi = $header;
    $gi =~ s/_.+//;
    
    my $taxid = $header;
    $taxid =~ s/.+_//;
    $taxid =~ s/;.+//; # some entries have more than 1 taxid - need to deal with this better at some point.... 

    ### use the subject sequences given by bsPrimerBlast to trim primer sequence aligned portion off sequence 

    ###
    ############  See Haikel's Nege2 primers, but the returned sequence/primers may need to be rev comp'd to find the match for trimming
    ### Negevirus_group2_F GTTGCWGGTCACGGTAARAC	Negevirus_group2_R CRTCAGCWGGAATWCGATAC

    my @seqsToTrimOff = @{ $seqTrimHash{$gi} };

    $seq = trimPrimers($seq, $seqsToTrimOff[0], $seqsToTrimOff[1], $gi);

    ### Get taxa info from database
    $sth->execute($taxid); 
    
    my ($species, $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $tax_name);
    while(my $row = $sth->fetchrow_hashref){
		$species      = "$row->{species_level}";
		$genus        = "$row->{genus_level}";
		$family       = "$row->{family_level}";
		$order        = "$row->{order_level}";
		$class        = "$row->{class_level}";
		$phylum       = "$row->{phylum_level}";
		$kingdom      = "$row->{kingdom_level}";
		$superkingdom = "$row->{superkingdom_level}";
		$tax_name     = "$row->{tax_name}";
    }


    ### print out 
    if(defined($species)) {
	# some of the entries in the tax db don't have species labeled >:-|
		if($species eq "NA") {
			$species = $tax_name;
		}

		$species =~ s/\'//g; # get rid of any ' in the species name - it messes up FastTree

		$speciesSeqCountHash{$species}++;
		
		# store taxa info for ampLen printing
		$taxaHash{$taxid} = join("\t", $superkingdom, $kingdom, $phylum, $class, $order, $family, $genus, $species);
		
		# store which taxid matches each gi for use when printing taxonomy table
		$giTaxaHash{$taxid . "\t" .
			    $gi . "\t" .
			    $taxaHash{$taxid}} = 1;
		
		# header with taxa info and unique number to make alignment program happy
		my $newHeader = ">"   . 
			$gi . ":" .
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
		if($speciesSeqCountHash{$species} <= $maxSeqPerSp) {
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

		$blastHitCount++; # store this for primer mismatch printing
	} else {
		print STDERR "Taxonomy not found for gi: $gi, taxid: $taxid\n";
		print STDERR "Make sure your blast database and taxonomy database " . 
		"were downloaded at about the same time\n";
	}
}

# Print out table with gi and taxa info for primerTree tree plotting
open my $giTaxaFile, ">", $outDir . "giTaxonomyFile.txt";
for my $giTaxaEntry (keys %giTaxaHash) {
	print $giTaxaFile $giTaxaEntry, "\n";
}
close $giTaxaFile;

# Print out information on how many mismatches there are
open my $mismatchFile, ">", $outDir . "primerMismatches.txt";
print $mismatchFile join("\t", "taxid", "superkingdom", "kingdom", "phylum", 
			 "class", "order", "family", "genus", "species", "direction", 
			 "mismatchTotal", "mismatch5Prime", "count", "totalCount"), "\n";
			 
open my $mismatchLocsFile, ">", $outDir . "primerMismatchLocs.txt";
print $mismatchLocsFile "# 1 is the 3 prime end of the primer\n";
print $mismatchLocsFile join("\t", "taxid", "superkingdom", "kingdom", "phylum", 
			 "class", "order", "family", "genus", "species", "direction", 
			 "mismatchLoc", "mismatchBase", "count", "totalCount"), "\n";

for my $taxid (keys %taxaHash) {
	# Mismatch counts				
	for my $mismatch (keys %{ $mismatchHash{$taxid} } ) {
	    print $mismatchFile $taxid . "\t" .
				$taxaHash{$taxid} . "\t" .
				$mismatch . "\t" .
				$mismatchHash{$taxid}{$mismatch} . "\t" .
				$blastHitCount, "\n"; 
	}
	# Mismatch locations	
	for my $mismatchLocs (keys %{ $mismatchLocsHash{$taxid} } ) {
	    my ($direction, $mismatchNums) = split("\t", $mismatchLocs); 

	    my @mismatchLocArray = split ",", $mismatchNums;
	    for my $loc (@mismatchLocArray) {
		my ($locNum, $locNt) = split "_", $loc;
		print $mismatchLocsFile $taxid . "\t" .
		    $taxaHash{$taxid} . "\t" .
		    $direction . "\t" .
		    $locNum . "\t" .
		    $locNt . "\t" .
		    $mismatchLocsHash{$taxid}{$mismatchLocs} . "\t" .
		    $blastHitCount, "\n"; 
	    }
	     
	    # For primers with no mismatches
	    if(scalar(@mismatchLocArray == 0)) {
		print $mismatchLocsFile $taxid . "\t" .
		    $taxaHash{$taxid} . "\t" .
		    $direction . "\t" .
		    "NA\tNA\t" .
		    $mismatchLocsHash{$taxid}{$mismatchLocs} . "\t" .
		    $blastHitCount, "\n"; 
	    }
	}
}

close $mismatchFile;
close $blastDbCmdResponse;
$dbh->disconnect;

##########################
### Print out amplicon length info with taxa
open my $ampliconLenFH, ">", $outDir . "ampliconLengths.txt";
print $ampliconLenFH join("\t", "taxid", "length", "count", "superkingdom", "kingdom", "phylum", 
			 "class", "order", "family", "genus", "species"), "\n";

for my $taxid (keys %ampLenHash) {
    if(exists($taxaHash{$taxid})) {
	for my $length (keys %{ $ampLenHash{$taxid} }) {
	    print $ampliconLenFH join("\t", $taxid, $length, $ampLenHash{$taxid}{$length}, $taxaHash{$taxid}), "\n";
	}
    } else {
	print STDERR "Taxa info for taxid ", $taxid, " not available during amplicon length file output\n";
    }
}

close $ampliconLenFH;

### Don't want to print out too many sequences for the alignment, so print out randomly if there are too many
if($maxAlignedSeqs <= scalar(@printable)) {
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

        ##### Need to keep in mind that the output fasta may  have multiple ">"s in the header - I don't think this is still true
$/ = "\n";



########################
### Print out a log of the number of taxa within each taxonomic level
open my $taxaSummaryFile, ">", $outDir . "taxaCountSummary.txt";
print $taxaSummaryFile "superkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tCount\n";

#for my $level ("subkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species") {
#    for my $taxName (keys %{ $taxNamesHash{$level} } ) {
#	print $taxaSummaryFile $level, "\t", $taxName, "\t", $taxNamesHash{$level}{$taxName}, "\n";
for my $taxName (keys %taxCountHash) {
    print $taxaSummaryFile $taxName, "\t", $taxCountHash{$taxName}, "\n";
    
}


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
### Calculate nucleotide distance between species within various levels
local $/ = "\n>";

open my $alignedFastaFH, $outDir . "seqsWithTaxaAligned.fasta" or die "Could not open aligned fasta output\n";
open my $distOutputFH, ">", $outDir . "distanceSummary.txt";

print $distOutputFH "CompLevel\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tMeanDist\tnCompared\n";

while(my $fasta = <$alignedFastaFH>) {
  chomp $fasta;
  my ($header,@sequences) = split "\n", $fasta;
  $header =~ s/>//;

  # remove these parts of the header to clean up the output
  #$header =~ s/^s-//;
  $header =~ s/:[kpcofgs]k*-/:/g;

  my $seq = join("", @sequences);
  $seq = lc($seq);

  my ($gi, $sp, $sk, $k, $p, $c, $o, $f, $g, undef) = split ":", $header;

  for my $fastaEntry (keys %alignedSeqHash) {
      # Only compare when either family, genus or species is the same between sequences
      if($alignedSeqHash{$fastaEntry}{family} eq $f || 
	     $alignedSeqHash{$fastaEntry}{genus} eq $g || 
	     $alignedSeqHash{$fastaEntry}{species} eq $sp) {

	  my $dist = dist($alignedSeqHash{$fastaEntry}{seq}, $seq);
	  # calc distance 
	  if($alignedSeqHash{$fastaEntry}{family} eq $f && $f ne "NA") {
	      my $keyToUse = "family\t" . 
			    $sk . "\t" . 
			    $k . "\t" . 
			    $p . "\t" . 
			    $c . "\t" . 
			    $o . "\t" . 
			    $f . "\t-\t-";

	      $distHash{$keyToUse}{sum} += $dist; 
	      # Check if $distHash entry exists, if not, initialize at 1
	      if(!exists($distHash{$keyToUse}{count})) {
		  $distHash{$keyToUse}{count} = 1;
	      }
	      $distHash{$keyToUse}{count}++;
	  }
	  if($alignedSeqHash{$fastaEntry}{genus} eq $g && $g ne "NA"){
	      my $keyToUse = "genus\t" . 
			    $sk . "\t" . 
			    $k . "\t" . 
			    $p . "\t" . 
			    $c . "\t" . 
			    $o . "\t" . 
			    $f . "\t" .
			    $g . "\t-";

	      $distHash{$keyToUse}{sum} += $dist; 
	      if(!exists($distHash{$keyToUse}{count})) {
		  $distHash{$keyToUse}{count} = 1;
	      }
	      $distHash{$keyToUse}{count}++;
	  }
	  if($alignedSeqHash{$fastaEntry}{species} eq $sp && $sp ne "NA"){
	      my $keyToUse = "species\t" . 
			    $sk . "\t" . 
			    $k . "\t" . 
			    $p . "\t" . 
			    $c . "\t" . 
			    $o . "\t" . 
			    $f . "\t" .
			    $g . "\t" . 
			    $sp;

	      $distHash{$keyToUse}{sum} += $dist; 
	      if(!exists($distHash{$keyToUse}{count})) {
		  $distHash{$keyToUse}{count} = 1;
	      }
	      $distHash{$keyToUse}{count}++;
	  }
      }
  }
  # store the new sequence for remaining comparisons
  $alignedSeqHash{$header}{seq} = $seq;
  $alignedSeqHash{$header}{family} = $f;
  $alignedSeqHash{$header}{genus} = $g;
  $alignedSeqHash{$header}{species} = $sp;
}

close $alignedFastaFH;

for my $distLevel (keys %distHash) {
    print $distOutputFH $distLevel, "\t", $distHash{$distLevel}{sum} / $distHash{$distLevel}{count}, "\t", $distHash{$distLevel}{count}, "\n";
}

close $distOutputFH;



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


sub trimPrimers {
    my $seq = $_[0];
    my $primer1 = $_[1];
    my $primer2 = $_[2];
    my $gi = $_[3];

    my $RCprimer1 = revComp($primer1);

    if($seq =~ /^$primer1/) {                               # primer 1 is at 5' end
	$seq =~ s/^$primer1//;
	$seq = trimPrimer2($seq, $primer2, $gi, "3prime");
    } elsif($seq =~ /^$RCprimer1/) {                        # primer 1 is at 5' end RC
	$seq =~ s/^$RCprimer1//;
	$seq = trimPrimer2($seq, $primer2, $gi, "3prime");
    } elsif($seq =~ /$primer1$/) {                          # primer 1 is at 3' end
	$seq =~ s/$primer1$//;
	$seq = trimPrimer2($seq, $primer2, $gi, "5prime");
    } elsif($seq =~ /$RCprimer1$/) {                        # primer 1 is at 3' end RC
	$seq =~ s/$RCprimer1$//;
	$seq = trimPrimer2($seq, $primer2, $gi, "5prime");
    } else {                                                # primer 1 not found
	print STDERR "First primer not found for trimming\n", $gi, "\t", $primer1, "\t", $seq, "\n\n";
	my $newSeq = trimPrimer2($seq, $primer2, $gi, "3prime");

	if(length($newSeq) != length($seq)) { # primer2 was trimmed
	    $seq = $newSeq;
	} else {                              # primer2 was not trimmed, try other end
	    $seq = trimPrimer2($seq, $primer2, $gi, "5prime");
	}
    }
    return $seq;
}

sub trimPrimer2 {
    my $seq = $_[0];
    my $primer2 = $_[1];
    my $gi = $_[2];
    my $end = $_[3];

    my $RCprimer2 = revComp($primer2);
	
    if($end eq "3prime") {
	if($seq =~ /$primer2$/) {
	    $seq =~ s/$primer2$//;
	} elsif($seq =~ /$RCprimer2$/) {
	    $seq =~ s/$RCprimer2$//;
	} else {
	    print STDERR "Second sequence not found for trimming\n\n";
	    print STDERR $gi, "\t", $primer2, "\t", $seq, "\n";
	}
    } else {
	if($seq =~ /^$primer2/) {
	    $seq =~ s/^$primer2//;
	} elsif($seq =~ /^$RCprimer2/) {
	    $seq =~ s/^$RCprimer2//;
	} else {
	    print STDERR "Second sequence not found for trimming\n\n";
	    print STDERR $gi, "\t", $primer2, "\t", $seq, "\n";
	}

    }
    return $seq;
}

sub dist {
    my $seq1 = $_[0];
    my $seq2 = $_[1];
    
    if(length($seq1) != length($seq2)) {
	print STDERR "compared sequences differ in length\n";
	die;
    }
    
    my $dist = 0;
    
    my @seq1Array = split("", $seq1);
    my @seq2Array = split("", $seq2);
    
    for(my $i = 0; $i < scalar(@seq1Array); $i++) {
	if($seq1Array[$i] ne $seq2Array[$i]) {
	    if($seq1Array[$i] !~ /[nN-]/ && $seq2Array[$i] !~ /[nN-]/) {
		$dist++;
	    }
	}
    }
    return $dist;
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
