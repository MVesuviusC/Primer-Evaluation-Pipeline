#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#'
#' @return
#' @export
#'
#' @examples
primer_mismatch_count <- function(output_dir, target_taxa, target_level) {
  primer_mismatches <- read.delim(paste(output_dir,
                                       "/bsPrimerTreeOut/primerMismatches.txt",
                                       sep = ""),
                                 stringsAsFactors = FALSE,
                                 header = TRUE)

  # Get rid of hits with "banned" words
  primer_mismatches <- as.data.frame(primer_mismatches[grep(banned_words,
                                                            primer_mismatches$species,
                                                            perl = T,
                                                            invert = T),])
  # Get rid of hits with only genus (no space in name)
  primer_mismatches <- primer_mismatches[grep(" ", primer_mismatches$species),]

  # Make forward and reverse better
  primer_mismatches$direction <- gsub("for",
                                      "Forward",
                                      primer_mismatches$direction)
  primer_mismatches$direction <- gsub("rev",
                                      "Reverse",
                                      primer_mismatches$direction)

  # Designate on/off-target
  primer_mismatches$OnTarget <- "Off-target"
  primer_mismatches$OnTarget[grepl(target_taxa,
                                   primer_mismatches[[target_level]])
                             ] <- "On-target"

  return(primer_mismatches)
}

#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#'
#' @return
#' @export
#'
#' @examples
#' test <- primer_mismatch_loc(output_dir = "C:/Users/MattDellLaptop/Desktop/mam16SOut2/",
#'                             target_taxa = "Mammalia",
#'                             target_level = "class")
primer_mismatch_loc <- function(output_dir, target_taxa, target_level) {
  primer_mismatch_locs <- read.delim(paste(output_dir,
                                         "/bsPrimerTreeOut/primerMismatchLocs.txt",
                                         sep = ""),
                                   stringsAsFactors = FALSE,
                                   comment.char = "#",
                                   header = TRUE)

  # Get rid of hits with "banned" words
  ####
  ### This makes the "totalCount" column incorrect since we're discarding some ***********************************************
  ### Need to correct this
  ####
  primer_mismatch_locs <- as.data.frame(primer_mismatch_locs[grep(banned_words,
                                                              primer_mismatch_locs$species,
                                                              perl = T,
                                                              invert = T),])

  # Get rid of hits with only genus (no space in name)
  primer_mismatch_locs <- primer_mismatch_locs[grep(" ", primer_mismatch_locs$species),]

  # Make forward and reverse better
  primer_mismatch_locs$direction <- gsub("for",
                                       "Forward",
                                       primer_mismatch_locs$direction)
  primer_mismatch_locs$direction <- gsub("rev",
                                       "Reverse",
                                       primer_mismatch_locs$direction)

  # Designate on/off-target
  primer_mismatch_locs$OnTarget <- "Off-target"
  primer_mismatch_locs$OnTarget[grepl(target_taxa,
                                      primer_mismatch_locs[[target_level]])] <- "On-target"

  # replace NA with 0 in
  return(primer_mismatch_locs)
}


#' Title
#'
#' @param output_dir
#' @param target_taxa
#' @param target_level
#' @param forward
#' @param reverse
#' @param target
#'
#' @return
#' @export
#'
#' @examples
#' test2 <- plot_primer_mismatch_locs(forward = "CGGTTGGGGTGACCTCGGA",
#'                                    reverse = "GCTGTTATCCCTAGGGTAACT",
#'                                    output_dir = "C:/Users/MattDellLaptop/Desktop/mam16SOut2/",
#'                                    target_taxa = "Mammalia",
#'                                    target_level = "class")
plot_primer_mismatch_locs <- function(forward, reverse, output_dir,
                                      target_taxa, target_level,
                                      target = "On-target") {

  primer_mismatch_locs <- primer_mismatch_loc(output_dir = output_dir,
                                               target_taxa = target_taxa,
                                               target_level = target_level) %>%
    dplyr::filter(OnTarget == target) %>%
    dplyr::filter(mismatchBase %in% c("A", "T", "G", "C"))

  # Make named list of primers to use as labels on figures
  Forward <- rev(strsplit(as.character(forward), split = "")[[1]])
  names(Forward) <- 1:length(Forward)

  Reverse <- rev(strsplit(as.character(reverse), split = "")[[1]])
  names(Reverse) <- 1:length(Reverse)

  plot_mismatch <- function(primer) {
    mismatch_plot <- primer_mismatch_locs %>%
        dplyr::filter(direction == primer) %>%
        dplyr::mutate(taxCount = length(unique(taxid))) %>%
        ggplot2::ggplot(., ggplot2::aes(x = mismatchLoc,
                                     y = count/taxCount,
                                     fill = mismatchBase)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = paste("Proportion of amplifiable targets with ",
                                    "primer mismatches at each location\n",
                                    target,
                                    " species, ",
                                    primer,
                                    " primer",
                                    sep = "")) +
        ggplot2::ylim(0,1) +
        ggplot2::xlab("5' end <-----------Primer position-----------> 3' end") +
        ggplot2::ylab("") +
        ggplot2::scale_x_reverse(limits = c(length(get(primer)), 1),
                                 breaks = 1:length(get(primer)),
                                 labels=get(primer))
  }

    # only create plot if there is data to use
    if(nrow(primer_mismatch_locs) > 0) {
      gridExtra::grid.arrange(grobs = lapply(c("Forward", "Reverse"),
                                             plot_mismatch))
    } else {
      warning("No data to plot for ", target, " ", primer)
    }
}
