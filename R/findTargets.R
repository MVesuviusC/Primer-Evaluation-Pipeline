#' Run initial blast searches
#'
#' @param forward forward primer
#' @param reverse reverse primer
#' @param blast_path path to blast version to use
#' @param blast_db path to blast database
#' @param tax_db path to taxonomy database created with make_tax_db()
#' @param threads number of threads to use
#' @param output_dir directory to write the files out to
#' @param assay_name name of the assay
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#' @param max_aligned_seqs Maximum number of sequences to keep for alignment
#'  and data processing
#'
#' @return none
#'
find_targets <- function(forward, reverse, assay_name, blast_path, blast_db,
                         tax_db, threads, output_dir, banned_words,
                         max_aligned_seqs) {

  package_dir <- find.package("bsPrimerTree")

  # Make output directory
  dir.create(path = output_dir,
             showWarnings = FALSE)

  # put together the command to run
  bsPrimerBlast_path <- paste(package_dir, "/exec/bsPrimerBlast.pl", sep = "")
  bsPrimerTree_path <- paste(package_dir, "/exec/bsPrimerTree.pl", sep = "")
  out_path <- paste(output_dir, "/bsPrimerTreeOut", sep = "")
  banned_words_with_quotes <- paste("\'", banned_words, "\'", sep = "")

  primer_blast_cmd <- paste("perl",
                            bsPrimerBlast_path,
                            "--tempDir", output_dir,
                            "--forward", forward,
                            "--reverse", reverse,
                            "--primerName", assay_name,
                            "--blastDb", blast_db,
                            "--blastVer", blast_path,
                            "--proc", threads,
                            "--minAmpLen 50",
                            "--maxAmpLen 1000",
                            "--verbose",
                            "|",
                            "perl",
                            bsPrimerTree_path,
                            "--inFile -",
                            "--blastDb", blast_db,
                            "--taxDb", tax_db,
                            "--outDir", out_path,
                            "--threads", threads,
                            "--maxSeqsPerSpecies 4",
                            "--maxAlignedSeqs", max_aligned_seqs,
                            "--verbose",
                            "--bannedWords", banned_words_with_quotes,
                            sep = " ")

  system(primer_blast_cmd)

  # check if the command succeeded

}

#' Load results from find_targets into output
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param output_dir character, directory to output intermediate files
#'
#' @return a list
#'
load_find_target_data <- function(bsPrimerTree, output_dir) {
  target_taxa <- bsPrimerTree$summary_table$target_taxa
  target_level <- bsPrimerTree$summary_table$target_level

  read_data <- function(in_file) {
    read.delim(paste(output_dir,
                     "/bsPrimerTreeOut/",
                     in_file,
                     sep = ""),
               stringsAsFactors = FALSE,
               comment.char = "#")
  }

  bsPrimerTree$amplicon_lengths <- amplicon_len(read_data("ampliconLengths.txt"),
                                                target_taxa = target_taxa,
                                                target_level = target_level)

  bsPrimerTree$distance_summary <- distance_data(read_data("distanceSummary.txt"),
                                                 target_taxa = target_taxa,
                                                 target_level = target_level)

  bsPrimerTree$acc_taxonomy <- read_data("accTaxonomyFile.txt")

  bsPrimerTree$primer_mismatches <-
    primer_mismatch_count(read_data("primerMismatches.txt"),
                          target_taxa = target_taxa,
                          target_level = target_level)

  bsPrimerTree$primer_mismatch_locs <-
    primer_mismatch_loc(read_data("primerMismatchLocs.txt"),
                        target_taxa = target_taxa,
                        target_level = target_level)

  bsPrimerTree$amplifiable <- read_data("taxaCountSummary.txt") %>%
    dplyr::select(-Count)

  bsPrimerTree$alignment <- ape::read.dna(paste(output_dir,
                                     "/bsPrimerTreeOut/seqsWithTaxaAligned.fasta",
                                     sep = ""), format = "fasta")

  bsPrimerTree$tree <- ape::read.tree(paste(output_dir,
                                "/bsPrimerTreeOut/tree.nwk",
                                sep = ""))

  bsPrimerTree
}

#' List on-target amplifiable species
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return a table of taxonomic data
#'
list_on_target_amplifiable <- function(bsPrimerTree, target_taxa, target_level) {
  # Get only on target hits
  on_target_amplifiable <- bsPrimerTree$amplifiable %>%
    dplyr::filter(get(target_level) == target_taxa)
}

#' List all known species within target taxonomic group
#'
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#' @param tax_db character, path to taxonomy database created with
#'   \code{\link{prepare_tax_db}}
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#'
#' @return A table of taxonomic data
#'
all_known_species <- function(target_taxa, target_level, tax_db, banned_words) {
  tax_name_with_quotes <- paste("\'", target_taxa, "\'", sep = "")
  tax_info <- paste(tax_name_with_quotes, target_level, sep = ",")

  known_taxa_cmd <-
    paste("getTaxaLocal.pl",
          "--dbName",  tax_db,
          "--taxName", tax_info,
          "--quiet")

  tax_data <- system(known_taxa_cmd, intern = TRUE)

  col_names <- strsplit(tax_data[1], split = "\t")

  potential_hit_df <- tax_data[2:length(tax_data)] %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = "value", sep = "\t", into = col_names[[1]], ) %>%
    dplyr::filter(grepl(banned_words, species) == FALSE)
}
