#' Title
#'
#' @param output_dir 
#' @param tax_db 
#'
#' @examples
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

#' Title
#'
#' @param output_dir 
#'
#'
#' @examples
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

#' Title
#'
#' @param output_dir 
#' @param target_taxa 
#' @param target_level 
#' @param forward 
#' @param reverse 
#'
#' @examples
potential_hits <- function(output_dir, target_taxa, target_level, 
                           forward, reverse, banned_words) {
  package_dir <- find.package("bsPrimerTree")
  
  # put together the command to run
  get_potential_path <- paste(package_dir, 
                         "/exec/getPotentialHits.pl", 
                         sep = "")
  input_file <- paste(output_dir, "/reblastResults.txt.gz", sep = "")
  taxa_file <- paste(output_dir, "/taxonomy.txt", sep = "")
  output_file <- paste(output_dir, "/reblastResults.txt.gz", sep = "")
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