#!/usr/bin/bash

makeblastdb \
  -in blastoSeqs.fasta \
  -parse_seqids \
  -blastdb_version 5 \
  -taxid_map test_map.txt \
  -title "test" \
  -dbtype nucl



# bsPrimerTree get seqs test
~/bin/ncbi-blast-2.10.1+/bin/blastdbcmd \
  -target_only \
  -db blastoSeqs.fasta \
  -entry_batch seqsToGet.txt \
  -outfmt ">%a_%T@%s" \
  > test.txt

# bsPrimerBlast test
blastn \
  -db blastoSeqs.fasta \
  -query blastoSeqs.fasta \
  -task blastn \
  -evalue 30000 \
  -word_size 7 \
  -max_hsps 100 \
  -max_target_seqs 50000 \
  -num_threads 2\
  -reward 1 \
  -penalty -1 \
  -gapopen 2 \
  -gapextend 1 \
  -outfmt "6 qseqid sgi qlen qstart qend sstart send slen sstrand qseq sseq staxids" \
  -ungapped \
  -sum_stats true \
  > test2.txt

# reblast test
blastn \
  -task blastn \
  -db blastoSeqs.fasta \
  -query blastoSeqs.fasta \
  -num_threads 2 \
  -outfmt "7 qseqid staxid score length qstart qend qlen sstart send slen sacc" \
  -max_hsps 1 \
  -max_target_seqs 10000 \
  > test3.txt

