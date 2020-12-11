#' Get genetic distance data on amplifiable sequences
#'
#' @param output_dir character, directory to output intermediate files
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#'
#' @return A table of taxa and distance data
#'
distance_data <- function(output_dir, target_taxa, target_level) {
  distance <- read.delim(paste(output_dir,
                               "/bsPrimerTreeOut/distanceSummary.txt",
                               sep = ""))

  distance$OnTarget <- "Off-target"
  distance$OnTarget[grepl(target_taxa, distance[[target_level]])] <- "On-target"

  distance_summary <- distance %>%
    dplyr::group_by(., CompLevel, OnTarget) %>%
    dplyr::mutate(., LevelAverage =
                    sum(MeanDist * nCompared) / sum(nCompared)) %>%
    dplyr::ungroup()

  return(distance_summary)
}

#' Plot genetic distance data between amplifiable sequences
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param levels_to_use character, taxonomic levels to plot - needs
#'   to be one of "skpcofgs"
#' @param target character, one of either On-target or Off-target
#'
#' @return
#' @export
#'
#' @examples
plot_distance <- function(bsPrimerTree,
                          levels_to_use = c("family", "genus", "species"),
                          target = "On-target") {
  output_dir <- bsPrimerTree$summary_table$output_dir
  target_taxa <- bsPrimerTree$summary_table$target_taxa
  target_level <- bsPrimerTree$summary_table$target_level

  distance_summary <- distance_data(output_dir = output_dir,
                                    target_taxa = target_taxa,
                                    target_level = target_level) %>%
    dplyr::filter(OnTarget == target)


  plot_dist <- function(level_in_use) {
    distance_summary %>%
      dplyr::filter(CompLevel == level_in_use) %>%
      ggplot2::ggplot(., ggplot2::aes(x = MeanDist,
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
