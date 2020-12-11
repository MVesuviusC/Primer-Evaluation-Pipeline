#' Get taxonomic data for blast results
#'
#' @param output_dir character, directory to output intermediate files
#' @param tax_db character, path to taxonomy database created with
#'   \code{\link{prepare_tax_db}}
#'
#' @return None, writes to a file
#'
get_taxa <- function(output_dir, tax_db) {
  package_dir <- find.package("bsPrimerTree")
  input_file <- paste(output_dir, "/reblastResults.txt.gz", sep = "")
  output_file <- paste(output_dir, "/taxonomy.txt", sep = "")
  prep_path <- paste(package_dir, "/exec/prepTaxids.pl", sep = "")

  get_taxa_cmd <-
    paste("perl",
          prep_path,
          "-i", input_file, "|",
          "getTaxaLocal.pl",
          "--dbName",  tax_db,
          "--taxids -",
          "--quiet",
          ">", output_file)
  system(get_taxa_cmd)
}

#' Quantify taxonomic specificity
#'
#' @param output_dir character, directory to output intermediate files
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#'
tax_specificity <- function(output_dir, banned_words) {
  package_dir <- find.package("bsPrimerTree")

  # put together the command to run
  find_top_path <- paste(package_dir,
                         "/exec/findTopBlastHitsWithTaxa.pl",
                         sep = "")
  count_num_path <- paste(package_dir, "/exec/countNumInRank.pl", sep = "")
  input_file <- paste(output_dir, "/reblastResults.txt.gz", sep = "")
  taxa_file <- paste(output_dir, "/taxonomy.txt", sep = "")
  banned_words_with_quotes <- paste("\'", banned_words, "\'", sep = "")

  find_top_cmd <- paste("perl", find_top_path,
                        "--blastIn", input_file,
                        "--taxa", taxa_file,
                        "--bannedWords", banned_words_with_quotes,
                        "|",
                        "perl", count_num_path,
                        "--input -")

  specificity_results <- system(find_top_cmd, intern = TRUE) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = "value", sep = "\t", into = c("label", "data")) %>%
    dplyr::pull(data, name = label)

  return(specificity_results)
}

#' Find sequences that could have been amplified
#'
#' @param output_dir character, directory to output intermediate files
#' @param forward character, forward primer
#' @param reverse character, reverse primer
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#'
potential_hits <- function(output_dir, target_taxa, target_level,
                           forward, reverse, banned_words) {
  package_dir <- find.package("bsPrimerTree")

  # put together the command to run
  get_potential_path <- paste(package_dir,
                         "/exec/getPotentialHits.pl",
                         sep = "")
  input_file <- paste(output_dir, "/reblastResults.txt.gz", sep = "")
  taxa_file <- paste(output_dir, "/taxonomy.txt", sep = "")
  banned_words_with_quotes <- paste("\'", banned_words, "\'", sep = "")

  get_potential_cmd <- paste("perl", get_potential_path,
                        "--alignVar 10",
                        "--blast", input_file,
                        "--taxa", taxa_file,
                        "--primerF", forward,
                        "--primerR", reverse,
                        "--bannedWords", banned_words_with_quotes)

  potential_hit_data <- system(get_potential_cmd, intern = TRUE)

  col_names <- strsplit(potential_hit_data[1], split = "\t")

  potential_hit_df <- potential_hit_data[2:length(potential_hit_data)] %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = "value", sep = "\t", into = col_names[[1]], ) %>%
    dplyr::filter(get(target_level) == target_taxa)

  return(potential_hit_df)
}
