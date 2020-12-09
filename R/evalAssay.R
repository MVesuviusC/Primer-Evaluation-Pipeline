#' Evaluate metabarcoding assay
#'
#' This function will take the provided information and evaluate the quality
#' of data that would be generated for a metabarcoding experiment.
#'
#' @details
#' Be aware that since this function uses blast, it may require a lot of
#' ram (40G+).
#'
#' @param forward character, forward primer
#' @param reverse character, reverse primer
#' @param blast_exe character, path to blast version to use
#' @param blast_db character, path to blast database
#' @param tax_db character, path to taxonomy database created with
#'   \code{\link{make_tax_db}}
#' @param output_dir character, directory to output intermediate files
#' @param threads numeric, number of processes to use
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#' @param assay_name character, name of the assay
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#' @param clean_up logical, if TRUE some of the unnecessary files will be
#'   removed at the end
#' @param ... other arguments
#'
#' @return A bsPrimerTree object, which is a list with the following elements,
#' \item{summary_table}{a list with summary metrics,
#'   most easily accessed through summary()}
#' \item{amplifiable}{A table of species that can be amplified by these primers}
#' \item{amplifiable_on_target}{A table of on-target species that can be
#'   amplified by these primers}
#' \item{species_seqd_at_locus}{A table of species that have sequence available
#'   in the blast nt database at this locus}
#' \item{missed_species}{A table of species that cannot be amplified by these
#'   primers}
#' \item{known_species}{A table of all known species within target group}
#' \item{species_unsequenced_at_locus}{A table of on-target species with no
#'   available sequences in the nt database}
#'
#' @export
#'
#' @examples
#'
#' test <- eval_assay(forward = "TGGTCGCAAGGCTGAAACTT",
#'                    reverse = "TTGCCTCCAGCTTCCCTACA",
#'                    output_dir = "testRun",
#'                    tax_db = "taxonomy.db",
#'                    assay_name = "testRun",
#'                    target_taxa = "Blastocystis",
#'                    target_level = "genus",
#'                    threads = 4)
#'
eval_assay <- function(forward, reverse, target_taxa, target_level, assay_name,
                       blast_exe = "blastn", blast_db = "nt", tax_db,
                       output_dir = "./output", threads = 1,
                       banned_words = banned_word_list, clean_up = FALSE, ...) {
  # Output
  output <- list()

  # Make sure target_taxa is Title Case and target_level is lowercase
  # and the taxa level is allowed
  if (!target_level %in% allowed_taxa_levels) {
    warning("The target taxonomic level (\"target_level\") provided in input
            options is not one of: ", paste(allowed_taxa_levels,
                                            collapse = ", ")
    )
    warning("Please modify target taxa and retry!")
    return()
  }

  # Set up output summary table
  output$summary_table <- list(assay_name = assay_name,
                               output_dir = output_dir,
                               primer_for = forward,
                               primer_rev = reverse,
                               target_taxa = target_taxa,
                               target_level = target_level,
                               banned_words = banned_words)

  depends_missing <- check_depends(blast_exe = blast_exe)

  if (depends_missing) {
    warning("Dependency missing. Stopping evaluation.", immediate. = TRUE)
    return()
  } else {
    # Run blast search on primer pairs
    find_targets(forward = forward,
                 reverse = reverse,
                 assay_name = assay_name,
                 blast_path = blast_exe,
                 blast_db = blast_db,
                 tax_db = tax_db,
                 threads = threads,
                 output_dir = output_dir,
                 banned_words = banned_words)

    # Add amplicon length info to table
    output$summary_table$median_on_target_amplicon_length <-
      amplicon_len(output_dir = output_dir,
                   target_taxa = target_taxa,
                   target_level = target_level) %>%
      dplyr::filter(onTarget == TRUE) %>%
      dplyr::pull(length) %>%
      median(na.rm = TRUE)

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

    # List all species amplifiable
    output$amplifiable <- list_all_amplifiable(output_dir = output_dir,
                                               target_taxa = target_taxa,
                                               target_level = target_level)

    # Count of all species that are amplifiable
    output$summary_table$speciesAmplifiableCount <-
      output$amplifiable %>%
      dplyr::pull(species) %>%
      unique() %>%
      length()

    # Get list of on-target amplifiable targets
    output$amplifiable_on_target <-
      list_on_target_amplifiable(output_dir = output_dir,
                                 target_taxa = target_taxa,
                                 target_level = target_level)

    # Count of on-target species that are amplifiable
    output$summary_table$onTargetSpeciesAmplifiableCount <-
      output$amplifiable_on_target %>%
      dplyr::pull(species) %>%
      unique() %>%
      length()

    # Percent of species amplifiable that are on-target
    output$summary_table$amplifiablePercentOnTarget <-
      output$summary_table$onTargetSpeciesAmplifiableCount /
      output$summary_table$speciesAmplifiableCount

    # Get info on distance between amplifiable targets
    ## Genus
    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachGenus <-
      distance_data(output_dir = output_dir,
                    target_taxa = target_taxa,
                    target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "genus") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)

    ## Species
    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachSpecies <-
      distance_data(output_dir = output_dir,
                    target_taxa = target_taxa,
                    target_level = target_level) %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "species") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)

    # Prep sequences for second blast
    warning("Starting blast on amplifiable sequences
            to get taxonmic specificity data.", immediate. = TRUE)
    reblast(blast_exe = blast_exe,
            blast_db = blast_db,
            threads = threads,
            output_dir = output_dir,
            target_taxa = target_taxa)

    # Get full taxonomy data for re-blast output
    get_taxa(output_dir = output_dir, tax_db = tax_db)

    # Quantify how many species/genera/families match each sequence queried
    # and what percentage of sequences match a single taxa
    output$summary_table <- c(output$summary_table,
                              tax_specificity(output_dir = output_dir,
                                              banned_words = banned_words))

    # Get list of potentially amplifiable species
    # that have locus sequence data available in NCBI
    output$species_seqd_at_locus <- potential_hits(output_dir = output_dir,
                                                   target_taxa = target_taxa,
                                                   target_level = target_level,
                                                   forward = forward,
                                                   reverse = reverse,
                                                   banned_words = banned_words)

    # Calculate the percent of species amplifiable
    output$summary_table$PercentAmplifiable <-
      (length(output$amplifiable_on_target$species) /
      length(output$species_seqd_at_locus$species)) * 100

    output$missed_species <- output$species_seqd_at_locus$species[
      output$species_seqd_at_locus$species %in%
      output$amplifiable_on_target$species]

    # Get list of all known taxa within target group
    output$known_species <-
      all_known_species(target_taxa = target_taxa,
                        target_level = target_level,
                        tax_db = tax_db,
                        banned_words = banned_words)

    # Calculate the percent of species with sequence for target locus
    output$summary_table$percentKnownSpeciesSeqd <-
      (length(output$species_seqd_at_locus$species) /
      length(output$species_in_database$species)) * 100

    # List missed species
    output$missed_species <-
      output$species_in_database[!output$species_in_database$species %in%
                                   output$amplifiable_on_target$species, ]

    # List unsequenced species
    output$species_unsequenced_at_locus <-
      output$species_in_database[!output$species_in_database$species %in%
                                   output$species_in_database$species, ]

    # Cleanup if I'm told to
    if (clean_up) {
      warning("Cleaning up some files", immediate. = TRUE)
      clean_up_cmd <- paste("rm ",
                            output_dir, "/reblastResults.txt.gz ",
                            output_dir, "/taxonomy.txt",
                            output_dir, "bsPrimerTreeOut/seqsWithTaxaAligned.fasta",
                            output_dir, "bsPrimerTreeOut/seqsWithTaxa.fasta",
                            sep = "")
      system(clean_up_cmd)
    }
  }
  class(output) <- "bsPrimerTree"
  return(output)
}

#' Plot all the plots
#'
#' One plot not enough? Make all the plots!
#'
#' @param bsPrimerTree a bsPrimerTree object returned by
#'   \code{\link{eval_assay}}
#'
#' @export
#'
#' @return Several plots
#'
#' @examples
#'
#' pdf("allThePlots.pdf", width = 10, height = 10)
#' plot_everything(blastoExample)
#' dev.off()
#'
plot_everything <- function(bsPrimerTree) {
  display_tree(bsPrimerTree = bsPrimerTree)
  display_wordcloud(bsPrimerTree = bsPrimerTree)
  plot_amplicon_len(bsPrimerTree = bsPrimerTree)
  plot_distance(bsPrimerTree = bsPrimerTree)
  plot_primer_mismatch_locs(bsPrimerTree = bsPrimerTree)
  plot_primer_mismatch_locs(bsPrimerTree = bsPrimerTree, target = "Off-target")
}

#' Check that the dependencies are present
#'
#' @param blast_exe character, path to blast executable
#'
#' @return logical, FALSE if all dependencies are present
#'
check_depends <- function(blast_exe = "blastn") {
  missing_dependency <- FALSE

  for (command_to_check in c(blast_exe, "blastdbcmd", "mafft", "perl",
                            "getTaxaLocal.pl", "FastTree")) {
    if (Sys.which(command_to_check) == "") {
      warning(command_to_check, " command not found. Please ensure this program
              is installed and in your $PATH")
      message(command_to_check, " command not found. Please ensure this program
              is installed and in your $PATH")
      missing_dependency <- TRUE
    }
  }

  # check if blastn has -sum_stats option
  check_blastn_cmd <- paste(blast_exe,
                            "-help",
                            "|",
                            "grep sum_stats")
  check_blastn <- suppressWarnings(system(check_blastn_cmd, intern = TRUE))

  if (length(check_blastn) == 0) {
    warning("Your version of blastn does not have the -sum_stats option.",
            immediate. = TRUE)
    warning("This version will not work for this pipeline.", immediate. = TRUE)
    warning("Either add a different blastn version to your $PATH or provide
            a path to the executable using the blast_exe option",
            immediate. = TRUE)

    missing_dependency <- TRUE
  }
  return(missing_dependency)
}

# This is a list of words that are used to exclude results
# with uncertain taxonomy
banned_word_list <- paste(c("\\ssp\\.",
                            "\\scf\\.",
                            "\\saff\\.",
                            "\\saffin\\.",
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
