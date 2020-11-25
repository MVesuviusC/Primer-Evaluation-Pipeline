#' Title
#'
#' @param output_dir 
#' @param tax_db 
#'
#' @return
#' @export
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

all_known_species <- function(target_taxa, target_level, tax_db) {
  
  #getTaxaLocal.pl --dbName ${taxonomyDb} --taxName ${targetTaxa},${targetLevel} > ${outDir}/couldHaveHit/knownSpecies.txt 
 
  package_dir <- find.package("bsPrimerTree")
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
    tibble::as.tibble() %>% 
    tidyr::separate(col = "value", sep = "\t", into = col_names[[1]], )
  
  return(potential_hit_df)
}