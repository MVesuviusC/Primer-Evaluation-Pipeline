#' Prepare taxonomy database to use with eval_assay
#'
#' This creates the taxonomy database used within bsPrimerTree.
#' This will download the data directly from NCBI and take roughly half an hour.
#' I recommend updating your blast database at the same time to avoid
#' inconsistencies between this database and the nt database.
#'
#' To use this, you need the code available at
#' https://github.com/MVesuviusC/getTaxa
#'
#' @param output_dir character, directory to output intermediate files
#' @param tax_db character, path to write taxonomy database
#'
#' @export
#'
#' @examples
#' \dontrun{
#' prepare_tax_db(output_dir = "\.", tax_db = "taxonomy.db")
#' }
prepare_tax_db <- function(output_dir, tax_db) {
  prep_cmd <-
    paste("makeTaxonomyDb.pl",
          "--outDir", output_dir,
          "--databaseName", tax_db,
          "--verbose")

  system(prep_cmd)
}
