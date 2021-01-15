#' Blast amplifiable sequences to assess taxonomic assignment
#'
#' @param output_dir character, directory to output intermediate files
#' @param blast_exe character, path to blast version to use
#' @param blast_db character, path to blast database
#' @param target_taxa character, taxonomic group targeted by the assay -
#'   needs to be one of "skpcofgs"
#' @param threads numeric, number of processes to use
#'
reblast <- function(output_dir, blast_exe, blast_db, target_taxa, threads) {
  unique_seqs(output_dir = output_dir, target_taxa = target_taxa)
  blast_files <- list.files(path = output_dir,
                            pattern = "seqsWithTaxaOnTargetUnique")

  blast <- function(file_name) {
    input_file <- paste(output_dir, "/", file_name, sep = "")
    output_file <- paste(output_dir, "/reblastResults.txt", sep = "")
    blast_cmd <-
      paste(blast_exe,
            "-task blastn",
            "-db", blast_db,
            "-query", input_file,
            "-num_threads", threads,
            "-outfmt \"7 qseqid staxid score length qstart qend qlen sstart send slen saccver\"",
            "-max_hsps 1",
            "-max_target_seqs 10000",
            ">>", output_file)

    system(blast_cmd)
  }
  lapply(blast_files, blast)
  system(paste("gzip -f ", output_dir, "/reblastResults.txt", sep = ""))
}

#' Get unique sequences
#'
#' @param output_dir character, directory to output intermediate files
#' @param target_taxa character, taxonomic group targeted by the assay -
#'   needs to be one of "skpcofgs"
#'
unique_seqs <- function(output_dir, target_taxa) {
  package_dir <- find.package("bsPrimerTree")

  # put together the command to run
  # Paths to files and scripts
  input_path <- paste(output_dir,
                      "/seqsWithTaxa.fasta",
                      sep = "")
  filterFastaNs_path <- paste(package_dir,
                              "/exec/filterFastaNs.pl",
                              sep = "")
  uniqueFastaBySpeciesSeq_path <- paste(package_dir,
                                        "/exec/uniqueFastaBySpeciesSeq.pl",
                                        sep = "")
  out_file_path <- paste(output_dir, "/seqsWithTaxaOnTargetUnique", sep = "")

  prep_seqs_cmd <- paste("perl",
                         filterFastaNs_path,
                         "--fasta", input_path,
                         "--targetTaxa", target_taxa,
                         "|",
                         "perl",
                         uniqueFastaBySpeciesSeq_path,
                         "--fasta -",
                         "|",
                         "split -a 3 -dl 1000 -",
                         out_file_path,
                         sep = " ")

  system(prep_seqs_cmd)
}
