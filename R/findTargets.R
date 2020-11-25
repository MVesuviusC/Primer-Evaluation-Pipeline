#' Title
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
#' @return
#'
#' @examples
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


#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#'
#' @return
#' @export
#'
#' @examples
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


