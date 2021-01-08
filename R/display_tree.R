#' Plot the a phylogenetic tree showing amplifiable targets
#' from data created from find_targets
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#' @param ... additional arguments passed to \code{\link{plot_tree_ranks}}
#'
#' @return a plot
#' @export
#'
#' @examples
#' \dontrun{
#' display_tree(bsPrimerTree = blasto_example)
#' }
display_tree <- function(bsPrimerTree, ...) {
  # need to trim the taxonomy data off of the tree labels
  # this data is kept because it is used in other parts of the analysis
  bsPrimerTree$tree$tip.label <- gsub(":.+", "", bsPrimerTree$tree$tip.label)

  print(primerTree::plot_tree_ranks(tree = bsPrimerTree$tree,
                                    taxonomy = bsPrimerTree$acc_taxonomy,
                                    ...))
}

#' Display a wordcloud of the taxonomic groups amplifiable
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#'
#' @return a plot
#' @export
#'
#' @examples
#' \dontrun{
#' display_wordcloud(bsPrimerTree = blasto_example)
#' }
display_wordcloud <- function(bsPrimerTree) {
  target_taxa <- bsPrimerTree$summary_table$target_taxa
  target_level <- bsPrimerTree$summary_table$target_level

  taxa_levels <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species")
  for_cloud <- bsPrimerTree$amplifiable %>%
    dplyr::mutate(Count = 1,
                  onTarget = dplyr::if_else(get(target_level) == !! target_taxa,
                                            "On-target",
                                            "Off-target",
                                            "Off-target")) %>%
    dplyr::select(-species) %>% # this adds too many entries
    tidyr::gather("Level", "Taxa", -Count, -onTarget) %>%
    dplyr::group_by(Taxa, Level, onTarget) %>%
    dplyr::summarize(Count = sum(Count)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Level) %>%
    dplyr::mutate(Percent = 100 * round(Count / sum(Count), 4)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Percent > 0.5) %>%
    dplyr::mutate(Level = factor(Level,
                                    levels = taxa_levels)) %>%
    dplyr::mutate(Taxa = tidyr::replace_na(Taxa, "ND"))

    print(ggplot2::ggplot(for_cloud,
                          ggplot2::aes(label = paste(Taxa,
                                                     " (", Percent, "%)",
                                                     sep = ""),
                                       size = Count,
                                       color = onTarget)) +
      ggwordcloud::geom_text_wordcloud_area(shape = "square",
                                            show.legend = TRUE,
                                            area_corr = FALSE) +
      ggplot2::scale_size_area(max_size = 5) +
      ggplot2::ggtitle("Taxa greater than 0.5%") +
      ggplot2::facet_wrap(~ Level, scales = "free") +
      ggplot2::guides(size = FALSE))
}
