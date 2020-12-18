#' Get genetic distance data on amplifiable sequences
#'
#' @param distance Distance table from bsPrimerTree.pl
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return A summary table of taxa and distance data
#'
distance_data <- function(distance, target_taxa, target_level) {

  distance$OnTarget <- "Off-target"
  distance$OnTarget[grepl(target_taxa, distance[[target_level]])] <- "On-target"

  distance_summary <- distance %>%
    dplyr::group_by(CompLevel, OnTarget) %>%
    dplyr::mutate(LevelAverage =
                    sum(MeanDist * nCompared) / sum(nCompared)) %>%
    dplyr::ungroup()
}

#' Plot genetic distance data between amplifiable sequences
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param levels_to_use character, taxonomic levels to plot - needs
#'   to be one of "skpcofgs"
#' @param target character, one of either On-target or Off-target
#'
#' @return a plot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_distance(blasto_example)
#' }
plot_distance <- function(bsPrimerTree,
                          levels_to_use = c("family", "genus", "species"),
                          target = "On-target") {
  target_taxa <- bsPrimerTree$summary_table$target_taxa
  target_level <- bsPrimerTree$summary_table$target_level

  plot_dist <- function(level_in_use) {
    bsPrimerTree$distance_summary %>%
      dplyr::filter(OnTarget == target) %>%
      dplyr::filter(CompLevel == level_in_use) %>%
      ggplot2::ggplot(ggplot2::aes(x = MeanDist,
                                      fill = OnTarget)) +
      ggplot2::geom_histogram(binwidth = 1) +
      ggplot2::ggtitle(paste(stringr::str_to_title(level_in_use),
                             "\nDistance between individual sequences\n",
                             target)) +
      ggplot2::xlab("") +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5))
  }

  gridExtra::grid.arrange(grobs = lapply(levels_to_use, plot_dist))
}
