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
  output_dir <- bsPrimerTree$summary_table$output_dir

  # Get taxa info
  gi_tax_df <- read.delim(paste(output_dir,
                                "/bsPrimerTreeOut/giTaxonomyFile.txt",
                                sep = ""),
                          header = TRUE)

  tree_data <- ape::read.tree(paste(output_dir,
                               "/bsPrimerTreeOut/tree.nwk",
                               sep = ""))

  # need to trim the taxonomy data off of the tree labels
  # this data is kept because it is used in other parts of the analysis
  tree_data$tip.label <- gsub(":.+", "", tree_data$tip.label)

  print(primerTree::plot_tree_ranks(tree = tree_data,
                                    taxonomy = gi_tax_df,
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
  output_dir <- bsPrimerTree$summary_table$output_dir
  target_taxa <- bsPrimerTree$summary_table$target_taxa
  target_level <- bsPrimerTree$summary_table$target_level

  # Get taxa info
  tax_summary <- read.delim(paste(output_dir,
                                "/bsPrimerTreeOut/taxaCountSummary.txt",
                                sep = ""),
                          header = TRUE,
                          stringsAsFactors = FALSE)

  for_cloud <- tax_summary %>%
    dplyr::mutate(., Count = 1,
                  onTarget = dplyr::if_else(get(target_level) == !! target_taxa,
                                            "On-target",
                                            "Off-target",
                                            "Off-target")) %>%
    dplyr::select(., -species) %>% # this adds too many entries
    tidyr::gather(., "Level", "Taxa", -Count, -onTarget) %>%
    dplyr::group_by(., Taxa, Level, onTarget) %>%
    dplyr::summarize(., Count = sum(Count)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(., Level) %>%
    dplyr::mutate(., Percent = 100 * round(Count / sum(Count), 4)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(., Percent > 0.5) %>%
    dplyr::mutate(., Level = factor(Level,
                                    levels = colnames(tax_summary[1:8]))) %>%
    dplyr::mutate(., Taxa = tidyr::replace_na(Taxa, "ND"))

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
