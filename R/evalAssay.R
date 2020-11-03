# output
output <- list()

#' Evaluate metabarcoding assay
#'
#' @param forward forward primer
#' @param reverse reverse primer
#' @param blast_path path to blast version to use
#' @param blast_db path to blast database
#' @param tax_db path to taxonomy database created with make_tax_db()
#' @param output_dir directory to output intermediate files
#' @param threads number of processes to use
#' @param target_taxa taxa targeted by the assay - needs to be one of "skpcofgs"
#' @param target_level taxonomic level of the targeted taxa
#' @param assay_name name of the assay
#' @param ... other arguments
#'
#' @return
#' @export
#'
#' @examples
#'
eval_assay <- function(forward, reverse, target_taxa, target_level, assay_name,
                       blast_path = "blastn", blast_db = "nt", tax_db,
                       output_dir = "./output/", threads = 1, ...) {

    # Make sure sure primers provided
  if (missing(forward) || missing(reverse)) {
    warning("Forward and reverse primers must be provided")
    return()
  }
  # Make sure taxa info provided
  if (missing(target_taxa) || missing(target_level)) {
    warning("Target taxa and taxa level must be provided")
    return()
  }

  # Make sure target_taxa is Title Case and target_level is lowercase
  # and the taxa level is allowed
  if(!target_level %in% allowed_taxa_levels) {
    warning("The target taxonomic level (\"target_level\") provided in input
            options is not one of: ", paste(allowed_taxa_levels,
                                            collapse = ", ")
            )

    warning("Please modify target taxa and retry!")
    return()
  }

  output$summary_table <- list(assay_name = assay_name,
                               output_dir = output_dir,
                               primer_for = forward,
                               primer_rev = reverse,
                               target_taxa = target_taxa,
                               target_level = target_level)

  depends_good <- check_depends()

  if (depends_good) {
    # Make directories for output
    dir.create(path = paste(output_dir,
                            "didHit",
                            sep = "/"),
               showWarnings = F,
               recursive = T)

    dir.create(path = paste(output_dir,
                            "reBlastOut",
                            sep = "/"),
               showWarnings = F)

    dir.create(path = paste(output_dir,
                            "couldHaveHit",
                            sep = "/"),
               showWarnings = F)

    find_targets(forward = forward,
                 reverse = reverse,
                 blast_path = blast_path,
                 blast_db = blast_db,
                 tax_db = tax_db,
                 threads = threads,
                 output_dir = output_dir)

    # Add amplicon length info to table
    output$summary_table$median_on-target_amplicon_length <-
      amplicon_len(output_dir = output_dir,
                   target_taxa = target_taxa,
                   target_level = target_level) %>%
      dplyr::filter(onTarget == TRUE) %>%
      dplyr::pull(length) %>%
      median(na.rm = TRUE)

    output$amplifiable <- list_on_target_amplifiable(output_dir = output_dir,
                                                     target_taxa = target_taxa,
                                                     target_level = target_level)

    # Get primer mismatch info
    output$summary_table$MeanOnTarget5PrimeMismatches <-
      primer_mismatch_count(output_dir = output_dir,
                            target_taxa = target_taxa,
                            target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target") %>%
      dplyr::pull(mismatch5Prime) %>%
      mean()

    output$summary_table$MeanOnTargetTotalMismatches <-
      primer_mismatch_count(output_dir = output_dir,
                            target_taxa = target_taxa,
                            target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target") %>%
      dplyr::pull(mismatchTotal) %>%
      mean()

    # Get info on distance between amplifiable targets
    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachGenus <-
      distance_data(output_dir = output_dir,
                    target_taxa = target_taxa,
                    target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "genus") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)

    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachSpecies <-
      distance_data(output_dir = output_dir,
                    target_taxa = target_taxa,
                    target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "species") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)






  } else {
    return()
  }
  class(output) <- "bsPrimerTree"
}


#' Check that the dependencies are present
#'
#' @return
#' @export
#'
#' @examples
#'
check_depends <- function() {


  return(TRUE)
}

# This is a list of words that are used to exclude results with uncertain taxonomy
banned_words <- paste(c(' sp\\.',
                        ' cf\\.',
                        ' aff\\.',
                        ' affin\\.',
                        "isolate",
                        "uncultured",
                        "symbiont",
                        "unidentified",
                        "unclassified",
                        "environmental"),
                      collapse = "|")

# List of taxonomic levels allowed
allowed_taxa_levels <- c("superkingdom",
                         "kingdom",
                         "phylum",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "species")
