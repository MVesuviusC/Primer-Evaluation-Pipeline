#!/bin/bash

# This script will download the NCBI blast databases needed as well as
# create the taxonomy database required
mkdir /usr/local/databases/blast

cd /usr/local/databases/blast/

update_blastdb.pl --decompress taxdb

update_blastdb.pl --decompress nt

cd /home/runDir

makeTaxonomyDb.pl -v --outDir /usr/local/databases/tax/ 
