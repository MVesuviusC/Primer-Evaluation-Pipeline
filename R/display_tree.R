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
