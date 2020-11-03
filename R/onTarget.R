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
on_target_amplifiable <- function(output_dir, target_taxa, target_level) {
  bsPrimerTree_hits <- read.delim(paste(output_dir,
                                        "/bsPrimerTreeOut/taxaCountSummary.txt",
                                        sep = ""),
                                  stringsAsFactors = FALSE,
                                  header = TRUE)

  # Make sure the provided taxaLevel argument is valid
  if(!target_level %in% allowed_taxa_levels) {
    warning("The target taxonomic level (\"target_level\") provided in
            input options is not one of: ", allowed_taxa_levels)
    warning("Please modify target taxa and retry!")
    knit_exit()
  }
  # Get only on target hits
  bsPrimerTree_hits <- subset(bsPrimerTree_hits, get(target_level) == target_taxa)

  return(bsPrimerTree_hits)
  # # Get rid of hits with "banned" words
  # bsPrimerTree_hits <- as.data.frame(bsPrimerTree_hits[grep(paste(bannedWords, collapse = "|"), bsPrimerTree_hits$species, perl = T, invert = T),])
  # # Get rid of hits with only genus (no space in name)
  # bsPrimerTree_hits <- bsPrimerTree_hits[grep(" ", bsPrimerTree_hits$species),]

}
