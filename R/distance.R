
#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#'
#' @return
#'
#' @examples
distance_data <- function(output_dir, target_taxa, target_level) {
  distance <- read.delim(paste(output_dir,
                               "/bsPrimerTreeOut/distanceSummary.txt",
                               sep = ""))

  distance$OnTarget <- "Off-target"
  distance$OnTarget[grepl(target_taxa, distance[[target_level]])] <- "On-target"

  # Get rid of hits with "banned" words
  distance <- as.data.frame(distance[grep(paste(banned_words, collapse = "|"),
                                          distance$species,
                                          perl = T,
                                          invert = T),])

  # Get rid of hits with only genus (no space in name)
  # Need to keep -
  distance <- distance[grep("[ ^-]", distance$species, perl = TRUE),]

  distance_summary <- distance %>%
    dplyr::group_by(., CompLevel, OnTarget) %>%
    dplyr::mutate(., LevelAverage = sum(MeanDist * nCompared) / sum(nCompared)) %>%
    dplyr::ungroup()

  return(distance_summary)
}


#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#' @param levels_to_use
#' @param target
#'
#' @return
#' @export
#'
#' @examples
plot_distance <- function(output_dir, target_taxa, target_level,
                          levels_to_use = c("family", "genus", "species"),
                          target = "On-target") {

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
