#' Plot the a phylogenetic tree showing amplifiable targets
#' from data created from find_targets
#'
#' @param output_dir directory where \code{\link{find_targets}} wrote output
#' @param ... additional arguments passed to \code{\link{plot_tree_ranks}}
#'
#' @return
#' @export
#'
#' @examples
display_tree <- function(output_dir, ...) {
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
  treeData$tip.label <- gsub(":.+", "", treeData$tip.label)

  primerTree::plot_tree_ranks(tree = tree_data, taxonomy = gi_tax_df, ...)
}


display_wordcloud <- function(output_dir, target_taxa, target_level) {
  # Get taxa info
  tax_summary <- read.delim(paste(output_dir,
                                "/bsPrimerTreeOut/taxaCountSummary.txt",
                                sep = "",),
                          header = TRUE,
                          stringsAsFactors = FALSE) 
  
  test <- tax_summary %>%
    dplyr::mutate(., Count = 1, 
                  onTarget = dplyr::if_else(get(target_level) == !! target_taxa, 
                                            "On-target", 
                                            "Off-target", 
                                            "Off-target")) %>%
    dplyr::select(., -species) %>% # this adds too many entries and they should all have Count = 1.... 
    tidyr::gather(., "Level", "Taxa", -Count, -onTarget) %>% 
    dplyr::group_by(., Taxa, Level, onTarget) %>% 
    dplyr::summarize(., Count = sum(Count)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(., Level) %>% 
    dplyr::mutate(., Percent = 100 * round(Count / sum(Count), 4)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(., Percent > 0.5) %>%
    dplyr::mutate(., Level = factor(Level, levels = colnames(tax_summary[1:8]))) %>%
    dplyr::mutate(., Taxa = tidyr::replace_na(Taxa, "ND"))
  
  ggplot(forCloud, aes(label = paste(Taxa, " (", Percent, "%)", sep = ""), 
                       size = Count, 
                       color = onTarget)) +
    geom_text_wordcloud_area(shape = "square", show.legend = TRUE) +
    scale_size_area(max_size = 3) +
    ggtitle("Taxa greater than 0.5%") +
    facet_wrap(~ Level, scales = "free") +
    guides(size = FALSE)  
}
