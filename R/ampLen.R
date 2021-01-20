#' Get data on the distribution of amplicon lengths for the tested assay
#'
#' @param amplicon_data character, amplicon data returned by bsPrimerTree.pl
#' @param target_taxa character, taxonomic group targeted by the assay -
#'   needs to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return
#' A table with taxonomy and amplicon length data
#'
amplicon_len <- function(amplicon_data, target_taxa, target_level) {
  # Mark which hits are on-target
  amplicon_data$onTarget <- grepl(target_taxa,
                                  amplicon_data[[target_level]])

  amplicon_data
}

#' Plot the distribution of amplicon lengths for the tested assay
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_amplicon_len(bsPrimerTree = blasto_example)
#' }
plot_amplicon_len <- function(bsPrimerTree) {
  amplicon_df <- bsPrimerTree$amplicon_lengths
  # In the current table, the count column is the number of hits
  # Plotting using this over represents commonly sequenced species
  amplicon_df$count <- 1
  bin_size <- 50
  print(ggplot2::ggplot(amplicon_df,
                        ggplot2::aes(x = as.numeric(length),
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
          ggplot2::scale_x_continuous(
            limits = c(0, max(amplicon_df$length + bin_size - 1)),
            breaks = round(seq(0,
                               max(amplicon_df$length + bin_size - 1),
                               by = bin_size),
                           1))
  )
}
