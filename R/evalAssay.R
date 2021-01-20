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
#'   \code{\link{prepare_tax_db}}
#' @param output_dir character, directory to output intermediate files
#' @param threads numeric, number of processes to use
#' @param target_taxa character, taxonomic group targeted by the assay - needs
#'   to be one of "skpcofgs"
#' @param target_level character, taxonomic level of the targeted taxa
#' @param assay_name character, name of the assay
#' @param banned_words list, a list of words used to identify sequences
#'   with uncertain taxonomy
#' @param max_aligned_seqs Maximum number of sequences to keep for alignment
#'  and data processing
#' @param num_permutations Maximum number of primer variants to include for
#'  primers with ambiguous bases
#' @param min_amp_len minimum length of amplicon to be included in evaluation
#' @param max_amp_len maximum length of amplicon to be included in evaluation
#' @param clean_up logical, if TRUE, intermediate files will be removed
#' @param ... other arguments
#'
#' @return A bsPrimerTree object, which is a list with the following elements,
#' \item{summary_table}{a list with summary metrics,
#'   most easily accessed through summary()}
#' \item{amplifiable}{A table of species that can be amplified by these primers}
#' \item{AmplifiableOnTarget}{A table of on-target species that can be
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
#' \dontrun{
#' blasto_example <- eval_assay(forward = "TGGTCGCAAGGCTGAAACTT",
#'                              reverse = "TTGCCTCCAGCTTCCCTACA",
#'                              output_dir = "testRun",
#'                              tax_db = "taxonomy.db",
#'                              assay_name = "blasto_example",
#'                              target_taxa = "Blastocystis",
#'                              target_level = "genus",
#'                              max_amp_len = 2000,
#'                              threads = 4)
#' }
eval_assay <- function(forward, reverse, target_taxa, target_level, assay_name,
                       blast_exe = "blastn", blast_db = "nt", tax_db,
                       output_dir = paste(tempdir(),
                                          "/",
                                          random_alphanumeric(20),
                                          sep = ""),
                       threads = 1, banned_words = banned_word_list,
                       max_aligned_seqs = 5000, num_permutations = 500,
                       min_amp_len = 0, max_amp_len = 2000, clean_up = TRUE, ...) {
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
                               primer_for = forward,
                               primer_rev = reverse,
                               target_taxa = target_taxa,
                               target_level = target_level,
                               banned_words = banned_words)

  depends_missing <- missing_depends(blast_exe = blast_exe)

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
                 banned_words = banned_words,
                 max_aligned_seqs = max_aligned_seqs,
                 num_permutations = num_permutations,
                 min_amp_len = min_amp_len,
                 max_amp_len = max_amp_len)

    # Load the files written out by bsPrimerTree into the output object
    output <- load_find_target_data(bsPrimerTree = output,
                                    output_dir = output_dir)

    # Add amplicon length info to summary table
    output$summary_table$MedianOnTargetAmpliconLength <-
      output$amplicon_lengths %>%
      dplyr::filter(onTarget == TRUE) %>%
      dplyr::pull(length) %>%
      median(na.rm = TRUE)

    # Add primer mismatch info to summary table
    output$summary_table$MeanOnTarget5PrimeMismatches <-
      output$primer_mismatches %>%
      dplyr::filter(OnTarget == "On-target") %>%
      dplyr::pull(mismatch5Prime) %>%
      mean(na.rm = TRUE)

    output$summary_table$MeanOnTargetTotalMismatches <-
      output$primer_mismatches %>%
      dplyr::filter(OnTarget == "On-target") %>%
      dplyr::pull(mismatchTotal) %>%
      mean(na.rm = TRUE)

    # Count of all species that are amplifiable
    output$summary_table$AllSpeciesAmplifiableCount <-
      output$amplifiable %>%
      dplyr::pull(species) %>%
      unique() %>%
      length()

    # Get list of on-target amplifiable targets
    output$amplifiableOnTarget <-
      list_on_target_amplifiable(output,
                                 target_taxa = target_taxa,
                                 target_level = target_level)

    # Count of on-target species that are amplifiable
    output$summary_table$OnTargetSpeciesAmplifiableCount <-
      output$amplifiableOnTarget %>%
      dplyr::pull(species) %>%
      unique() %>%
      length()

    # Percent of species amplifiable that are on-target
    output$summary_table$AmplifiablePercentOnTarget <-
      (output$summary_table$OnTargetSpeciesAmplifiableCount /
      output$summary_table$AllSpeciesAmplifiableCount) * 100

    # Get info on distance between amplifiable targets
    ## Genus
    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachGenus <-
      output$distance_summary %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "genus") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)

    ## Species
    output$summary_table$MeanOnTargetDistBetweenSeqsWithinEachSpecies <-
      output$distance_summary %>%
      dplyr::filter(OnTarget == "On-target", CompLevel == "species") %>%
      dplyr::pull(LevelAverage) %>%
      head(n = 1)

    # Prep sequences for second blast
    message("Starting blast on amplifiable sequences
            to get taxonmic specificity data.")
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

    # Get list of all known taxa within target group
    output$known_species <-
      all_known_species(target_taxa = target_taxa,
                        target_level = target_level,
                        tax_db = tax_db,
                        banned_words = banned_words)

    # Calculate the percent of species with sequence for target locus
    output$summary_table$PercentKnownSpeciesSeqd <-
      length(output$species_seqd_at_locus$species) /
      length(output$known_species$species)

    # List missed species
    output$missed_species <-
      output$known_species[!output$known_species$species %in%
                                   output$amplifiableOnTarget$species, ]

    # List unsequenced species
    output$species_unsequenced_at_locus <-
      output$known_species[!output$known_species$species %in%
                                   output$species_seqd_at_locus$species, ]

    # Cleanup if I'm told to
    if (clean_up) {
      message("Cleaning up some files")
      clean_up_cmd <- paste("rm ",
                            output_dir, "/*",
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
#' \dontrun{
#' pdf("allThePlots.pdf", width = 10, height = 10)
#' plot_everything(blasto_example)
#' dev.off()
#' }
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
missing_depends <- function(blast_exe = "blastn") {
  missing_dependency <- FALSE

  for (command_to_check in c(blast_exe, "blastdbcmd", "mafft", "perl",
                            "getTaxaLocal.pl", "FastTree")) {
    if (Sys.which(command_to_check) == "") {
      message(command_to_check, " command not found. Please ensure this program
              is installed and in your $PATH")
      message(command_to_check, " command not found. Please ensure this program
              is installed and in your $PATH")
      missing_dependency <- TRUE
    }
  }

  # check if blastn has -sum_stats option
  if (Sys.which(blast_exe) != "") {
    check_blastn_cmd <- paste(blast_exe,
                              "-help",
                              "|",
                              "grep sum_stats")
    check_blastn <- suppressWarnings(system(check_blastn_cmd, intern = TRUE))

    if (length(check_blastn) == 0) {
      message("Your version of blastn does not have the -sum_stats option.",
              immediate. = TRUE)
      message("This version will not work for this pipeline.", immediate. = TRUE)
      message("Either add a different blastn version to your $PATH or provide
            a path to the executable using the blast_exe option",
              immediate. = TRUE)

      missing_dependency <- TRUE
    }
  }

  missing_dependency
}


#' Summarize bsPrimerTree object
#'
#' @param object A bsPrimerTree object
#' @param ... Ignored options
#'
#' @return A bsPrimerTree object from \code{\link{eval_assay}}
#' @export
#'
#' @examples
#' \dontrun{
#' summary(blasto_example)
#' }
summary.bsPrimerTree <- function(object, ...) {
  summary_data <- stack(object$summary_table)[, c(2, 1)]
  colnames(summary_data) <- c("Label", object$summary_table$assay_name)

  package_dir <- find.package("bsPrimerTree")
  field_descriptions <- read.delim(paste(package_dir,
                                         "/summaryTableFieldDescriptions.txt",
                                         sep = ""),
                                   header = TRUE,
                                   stringsAsFactors = FALSE)

  summary_output <- dplyr::full_join(field_descriptions,
                                     summary_data,
                                     by = "Label") %>%
    dplyr::relocate(Label) %>%
    dplyr::relocate(Description, .after = dplyr::last_col())

  summary_output
}

# This is a list of words that are used to exclude results
# with uncertain taxonomy
# I may make these into functions to get the data so I can make them private
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

#' bsPrimerTree results for the mammalian 16S primers
#' @name mammal_example
#' @docType data
NULL

#' PrimerTree results for the blastocystis 18S primers
#' @name blasto_example
#' @docType data
NULL
