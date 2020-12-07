#' Get data on the distribution of amplicon lengths for the tested assay
#'
#' @param output_dir directory to output intermediate files
#' @param target_taxa taxonomic group targeted by the assay - needs to be one of "skpcofgs"
#' @param target_level taxonomic level of the targeted taxa
#'
#' @return
#' a dataframe with taxonomy and amplicon length data
#' 
#' @export
#'
#' @examples
#' amp_lens <- amplicon_len(output_dir = "output",
#'                          target_taxa = "Blastocystis", 
#'                          target_level = "genus")
#'                          
amplicon_len <- function(output_dir, target_taxa, target_level) {
  # Pull in data
  amplicon_len_df <- read.delim(paste(output_dir,
                                      "/bsPrimerTreeOut/ampliconLengths.txt",
                                      sep = ""),
                                header = TRUE)

  # Mark which hits are on-target
  amplicon_len_df$onTarget <- grepl(target_taxa,
                                    amplicon_len_df[[target_level]])

  return(amplicon_len_df)
}


#' Plot the distribution of amplicon lengths for the tested assay
#'
#' @param output_dir directory to output intermediate files
#' @param target_taxa taxonomic group targeted by the assay - needs to be one of "skpcofgs"
#' @param target_level taxonomic level of the targeted taxa
#'
#' @export
#'
#' @examples
#' plot_amplicon_len(output_dir = "output",
#'                   target_taxa = "Blastocystis", 
#'                  target_level = "genus")
#'                  
plot_amplicon_len <- function(output_dir, target_taxa, target_level) {
  amplicon_df <- amplicon_len(output_dir, target_taxa, target_level)

  print(ggplot2::ggplot(amplicon_df, ggplot2::aes(x = as.numeric(length), 
                                                  y = count, 
                                                  fill = onTarget)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::ggtitle("Amplicon Lengths") +
          ggplot2::xlab("Amplicon Length (bp)") +
          ggplot2::ylab("Count") + 
          ggplot2::facet_wrap(~ onTarget, 
                              ncol = 1, 
                              scales = "free_y") +
          ggplot2::labs(fill = "On target") +
          ggplot2::scale_x_continuous(breaks = round(seq(0, 
                                                         max(amplicon_df$length), 
                                                         by = 50), 
                                                     1))
  )
}
