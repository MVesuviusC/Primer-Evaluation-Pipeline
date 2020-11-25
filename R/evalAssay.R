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
#' @param banned_words 
#'
#' @return
#' @export
#'
#' @examples
#'
eval_assay <- function(forward, reverse, target_taxa, target_level, assay_name,
                       blast_path = "blastn", blast_db = "nt", tax_db,
                       output_dir = "./output/", threads = 1, 
                       banned_words = banned_word_list, clean_up = FALSE, ...) {
  
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
                               target_level = target_level,
                               banned_words = banned_words)
  
  depends_missing <- check_depends(blast_path = blast_path)
  
  if (depends_missing) {
    warning("Dependency missing. Stopping evaluation.", immediate. = TRUE)
    return()
  } else {
    # Run blast search on primer pairs
    find_targets(forward = forward,
                 reverse = reverse,
                 assay_name = assay_name,
                 blast_path = blast_path,
                 blast_db = blast_db,
                 tax_db = tax_db,
                 threads = threads,
                 output_dir = output_dir,
                 banned_words = banned_words)
    
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
    
    # Get list of on-target amplifiable targets
    output$amplifiable_on_target <- list_on_target_amplifiable(output_dir = output_dir,
                                                          target_taxa = target_taxa,
                                                          target_level = target_level)
    
    output$summary_table$speciesAmplifiableCount <- output$amplifiable_on_target %>%
      dplyr::pull(species) %>%
      length(uniq(.))
    
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
    
    # Prep sequences for second blast
    reblast(blast_path = blast_path,
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
    # - have locus sequence data available in NCBI
    output$potentially_amplifiable <- potential_hits(output_dir = output_dir,
                                                     target_taxa= target_taxa, 
                                                     target_level= target_level, 
                                                     forward = forward, 
                                                     reverse = reverse,
                                                     banned_words = banned_words)
    
    # Calculate the percent of species amplifiable
    output$summary_table$PercentAmplifiable <- 
      (length(output$amplifiable_on_target$species) /
      length(output$potentially_amplifiable$species)) * 100
    
    output$missed_species <- output$potentially_amplifiable$species[
      output$potentially_amplifiable$species %in% 
      output$amplifiable_on_target$species]
    
    # Get list of all known taxa within target group
    outtput$summary_table$species_in_database <- 
      all_known_species(target_taxa = target_taxa, 
                        target_level = target_level, 
                        tax_db = tax_db)
    
    # Calculate the percent of species with sequence for target locus
    percent_known_species_seqd <- 
      (length(output$potentially_amplifiable$species) /
      length(output$species_in_database$species)) * 100
    
    
    # Cleanup if I'm told to
    if(clean_up) {
      clean_up_cmd <- paste("rm ",
                            output_dir, "/reblastResults.txt.gz", 
                            sep = "")
      system(clean_up_cmd)
    }
    
    # List missed species
    
    
    # List unsequenced species
    
  }
  class(output) <- "bsPrimerTree"
}


#' Check that the dependencies are present
#'
#' @param blast_path 
#'
#' @return
#' @export
#'
#' @examples
#'
check_depends <- function(blast_path) {
  # c("blastn", "blastdbcmd", "mafft", "perl", "getTaxa.pl")

  missing_dependency <- FALSE
  
  for(command_to_check in c(blast_path, "blastdbcmd", "mafft", "perl", 
                            "getTaxaLocal.pl")) {
    if(Sys.which(command_to_check) == "") {
      warning(command_to_check, " command not found. Please ensure this program 
              is installed and in your $PATH")
      message(command_to_check, " command not found. Please ensure this program 
              is installed and in your $PATH")
      missing_dependency <- TRUE
    }
  }
  
  # check if blastn has -sum_stats option
  check_blastn_cmd <- paste("blastn",
                            "-help",
                            "|",
                            "grep sum_stats")
  check_blastn <- suppressWarnings(system(check_blastn_cmd, intern = TRUE))
  
  if(length(check_blastn) == 0) {
    warning("Your version of blastn does not have the -sum_stats option.", 
            immediate. = TRUE)
    warning("This version will not work for this pipeline.", immediate. = TRUE)
    warning("Either add a different blastn version to your $PATH or provide 
            a path to the executable using the blast_path option", 
            immediate. = TRUE)
    
    missing_dependency = TRUE
  }
  return(missing_dependency)
}

# This is a list of words that are used to exclude results with uncertain taxonomy
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
