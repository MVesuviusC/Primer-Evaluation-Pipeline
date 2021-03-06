% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primerMismatch.R
\name{primer_mismatch_loc}
\alias{primer_mismatch_loc}
\title{Add info to location of primer nucleotide mismatches table}
\usage{
primer_mismatch_loc(mismatch_locs, target_taxa, target_level)
}
\arguments{
\item{mismatch_locs}{Primer mismatch location data from bsPrimerTree.pl}

\item{target_taxa}{character, taxonomic group targeted by the assay - needs
to be one of "skpcofgs"}

\item{target_level}{character, taxonomic level of the targeted taxa}
}
\value{
A table of primer mismatch location data
}
\description{
Add info to location of primer nucleotide mismatches table
}
