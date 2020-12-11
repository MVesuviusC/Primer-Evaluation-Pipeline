#' Run initial blast searches
#'
#' @param forward forward primer
#' @param reverse reverse primer
#' @param blast_path path to blast version to use
#' @param blast_db path to blast database
#' @param tax_db path to taxonomy database created with make_tax_db()
#' @param threads number of threads to use
#' @param output_dir directory to write the files out to
#' @param assay_name
#'
#' @return none
#'
find_targets <- function(forward, reverse, assay_name, blast_path, blast_db,
                         tax_db, threads, output_dir, banned_words) {

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
                            "--forward",  forward,
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
                            "--maxAlignedSeqs 5000",
                            "--verbose",
                            "--bannedWords", banned_words_with_quotes,
                            sep = " ")

  system(primer_blast_cmd)

  # check if the command succeeded

}

#' List on-target amplifiable species
#'
#' @param output_dir character, directory to output intermediate files
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return a table of taxonomic data
#'
list_on_target_amplifiable <- function(output_dir, target_taxa, target_level) {
  # get list of species amplified
  bsPrimerTree_hits <- read.delim(paste(output_dir,
                                        "/bsPrimerTreeOut/taxaCountSummary.txt",
                                        sep = ""))

  # Get only on target hits
  bsPrimerTree_hits <- bsPrimerTree_hits %>%
    dplyr::filter(get(target_level) == target_taxa) %>%
    dplyr::select(-Count)

  return(bsPrimerTree_hits)
}

#' List all amplifiable species
#'
#' @param output_dir character, directory to output intermediate files
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return A table of taxonomic data
#'
list_all_amplifiable <- function(output_dir, target_taxa, target_level) {
  # get list of species amplified
  bsPrimerTree_hits <- read.delim(paste(output_dir,
                                        "/bsPrimerTreeOut/taxaCountSummary.txt",
                                        sep = ""))

  # Get only on target hits
  bsPrimerTree_hits <- bsPrimerTree_hits %>%
    dplyr::select(-Count)

  return(bsPrimerTree_hits)
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
    dplyr::filter(., grepl(banned_words, species) == FALSE)

  return(potential_hit_df)
}
