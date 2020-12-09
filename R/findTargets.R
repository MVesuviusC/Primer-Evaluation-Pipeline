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
#'   \code{\link{make_tax_db}}
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
