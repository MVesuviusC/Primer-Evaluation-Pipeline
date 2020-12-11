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
  print(ggplot2::ggplot(bsPrimerTree$amplicon_lengths,
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
          ggplot2::scale_x_continuous(breaks = round(seq(0,
                                                         max(amplicon_df$length),
                                                         by = 50),
                                                     1))
  )
}
