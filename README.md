# Primer evaluation

## To do

   * Fix primer mismatch total issue with removing banned words
   * Kick out genus only systematically


## Overview



## Installation

### Dependencies

Note that this pipeline likely requires around 40G of RAM minimum to run

## Usage



## Output



## Docker container


## Other scripts

### bsPrimerBlast
	and
### bsPrimerTree

This code is very much a work in progress and will likely change a good bit over time. Use with this in mind. If you try it out, feel free to send feedback to matthewvc1@gmail.com

The documentation is so far pretty much non-existent but will improve once this code is closer to being ready for general use.

The point of this code is to make a locally available version of the primer specificity assessment portion of primerBlast. I've copied the algorithm, with a few minor changes, from the primerBlast paper and code. 


Dependencies:
 - Dendroscope
 - getTaxa.pl
 - FastTree
 - Blastn - not version 2.8.1
   - blast nt database
 - Mafft
 - Imagemagick - optional - to make jpg from svg

