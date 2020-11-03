

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
#'
amplicon_len <- function(output_dir, target_taxa, target_level) {
  # Pull in data
  amplicon_len_df <- read.delim(paste(output_dir,
                                      "/bsPrimerTreeOut/ampliconLengths.txt",
                                      sep = ""),
                                header = TRUE)

  # Get rid of hits with "banned" words
  ## find location of species with banned words
  bad_species_locs <- grep(banned_words,
                           amplicon_len_df$species,
                           perl = T,
                           invert = T)
  ## remove them
  amplicon_len_df <- as.data.frame(amplicon_len_df[bad_species_locs,])

  # Get rid of hits with only genus (no space in name)
  amplicon_len_df <- amplicon_len_df[grep(" ", amplicon_len_df$species),]

  # Mark which hits are on-target
  amplicon_len_df$onTarget <- grepl(target_taxa,
                                    amplicon_len_df[[target_level]])

  return(amplicon_len_df)
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
plot_amplicon_len <- function(output_dir, target_taxa, target_level) {
  amplicon_df <- amplicon_len(output_dir, target_taxa, target_level)

  print(ggplot(amplicon_df, aes(x = as.numeric(length), y = count, fill = onTarget)) +
    geom_bar(stat = "identity") +
    ggtitle("Amplicon Lengths") +
    xlab("Amplicon Length (bp)") +
    ylab("Count") + facet_wrap(~ onTarget, ncol = 1, scales = "free_y") +
    labs(fill = "On target") +
    scale_x_continuous(breaks = round(seq(0, max(amplicon_df$length), by = 50), 1))
  )
}
