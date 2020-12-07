#' Prepare taxonomy database to use with eval_assay
#' 
#' This creates the taxonomy database used within bsPrimerTree.
#' This will download the data directly from NCBI and take roughly half an hour.
#' I recommend updating your blast database at the same time to avoid 
#' inconsistencies between this database and the nt database.
#'
#' @param output_dir 
#' @param tax_db 
#'
#' @export
#'
#' @examples
#' 
#' prepare_tax_db(output_dir = "\.", tax_db = "taxonomy.db")
prepare_tax_db <- function(output_dir, tax_db) {
  prep_cmd <- 
    paste("makeTaxonomyDb.pl",
          "--outDir", output_dir,
          "--databaseName", tax_db,
          "--verbose")
  
  system(prep_cmd)
}