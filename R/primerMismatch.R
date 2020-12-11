#' Add info to primer nucleotide mismatches table
#'
#' @param mismatch_data mismatch data output by bsPrimerTree.pl
#' @param target_taxa
#' @param target_level
#'
#' @return A table of primer mismatch data
#'
primer_mismatch_count <- function(mismatch_data, target_taxa, target_level) {
  # Make forward and reverse better
  mismatch_data$direction <- gsub("for",
                                      "Forward",
                                      mismatch_data$direction)
  mismatch_data$direction <- gsub("rev",
                                      "Reverse",
                                      mismatch_data$direction)

  # Designate on/off-target
  mismatch_data$OnTarget <- "Off-target"
  mismatch_data$OnTarget[grepl(target_taxa,
                                   mismatch_data[[target_level]])
                             ] <- "On-target"

  mismatch_data
}

#' Add info to location of primer nucleotide mismatches table
#'
#' @param mismatch_locs Primer mismatch location data from bsPrimerTree.pl
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return A table of primer mismatch location data
#'
primer_mismatch_loc <- function(mismatch_locs, target_taxa, target_level) {
  # Make forward and reverse better
  mismatch_locs$direction <- gsub("for",
                                       "Forward",
                                       mismatch_locs$direction)
  mismatch_locs$direction <- gsub("rev",
                                       "Reverse",
                                       mismatch_locs$direction)

  # Designate on/off-target
  mismatch_locs$OnTarget <- "Off-target"
  mismatch_locs$OnTarget[grepl(target_taxa,
                                      mismatch_locs[[target_level]])] <-
    "On-target"

  mismatch_locs
}

#' Plot primer nucleotide mismatch locations
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param target either "On-target" or "Off-target"
#'
#' @return a plot of mismatch locations
#' @export
#'
#' @examples
#' \dontrun{
#' plot_primer_mismatch_locs(bsPrimerTree = blasto_example,
#'                                    target = "On-target")
#' }
plot_primer_mismatch_locs <- function(bsPrimerTree, target = "On-target") {
  forward <- bsPrimerTree$summary_table$primer_for
  reverse <- bsPrimerTree$summary_table$primer_rev

  # Make named list of primers to use as labels on figures
  Forward <- rev(strsplit(as.character(forward), split = "")[[1]])
  names(Forward) <- seq_len(length(Forward))

  Reverse <- rev(strsplit(as.character(reverse), split = "")[[1]])
  names(Reverse) <- seq_len(length(Reverse))

  plot_mismatch <- function(primer) {
    mismatch_plot <- bsPrimerTree$primer_mismatch_locs %>%
      dplyr::filter(OnTarget == target) %>%
      dplyr::filter(mismatchBase %in% c("A", "T", "G", "C", NA)) %>%
        dplyr::filter(direction == primer) %>%
        dplyr::mutate(taxCount = length(unique(taxid))) %>%
        ggplot2::ggplot(., ggplot2::aes(x = mismatchLoc,
                                     y = count / taxCount,
                                     fill = mismatchBase)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = paste("Proportion of amplifiable targets with ",
                                    "primer mismatches at each location\n",
                                    target,
                                    " species, ",
                                    primer,
                                    " primer",
                                    sep = "")) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("5' end <-----------Primer position-----------> 3' end") +
        ggplot2::ylab("") +
        ggplot2::scale_x_reverse(limits = c(length(get(primer)), 1),
                                 breaks = seq_len(length(get(primer))),
                                 labels = get(primer))
  }
  gridExtra::grid.arrange(grobs = lapply(c("Forward", "Reverse"),
                                         plot_mismatch))
}
