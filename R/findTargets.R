

#' Title
#'
#' @param forward forward primer
#' @param reverse reverse primer
#' @param blast_path path to blast version to use
#' @param blast_db path to blast database
#' @param tax_db path to taxonomy database created with make_tax_db()
#' @param threads number of threads to use
#' @param output_dir directory to write the files out to
#'
#' @return
#'
#' @examples
#'
find_targets <- function(forward, reverse, blast_path, blast_db,
                         tax_db, threads, output_dir) {

  package_dir <- find.package("bsPrimerTree")

  # put together the command to run
  primer_blast_cmd <- paste("perl ",
                            package_dir, "/exec/countNumInRank.pl ",
                            "--forward ",  forward,
                            "--reverse ", reverse,
                            "--blastDb ", blast_db,
                            "--blastVer ", blast_path,
                            "--proc ", threads,
                            "--minAmpLen 50 ",
                            "--maxAmpLen 1000 ",
                            "| ",
                            "perl ",
                            package_dir, "/exec/bsPrimerTree.pl ",
                            "--inFile - ",
                            "--blastDb ", blast_db,
                            "--taxDb ", tax_db,
                            "--outDir ", output_dir, "/bsPrimerTreeOut ",
                            "--threads ", threads,
                            "--maxSeqsPerSpecies 4 ",
                            "--maxAlignedSeqs 5000 ",
                            sep = "")

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

  # Get only data from on-target hits
  bsPrimerTree_hits <- subset(bsPrimerTree_hits, get(target_level) == target_taxa)

  # Get rid of hits with "banned" words
  bsPrimerTree_hits <- bsPrimerTree_hits[grep(banned_words,
                                              bsPrimerTree_hits$species,
                                              perl = T,
                                              invert = T),]

  # Get rid of hits with only genus (no space in name)
  bsPrimerTree_hits <- bsPrimerTree_hits[grep(" ", bsPrimerTree_hits$species),]

  return(bsPrimerTree_hits)
}


