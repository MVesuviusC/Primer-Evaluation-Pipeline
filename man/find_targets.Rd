% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findTargets.R
\name{find_targets}
\alias{find_targets}
\title{Run initial blast searches}
\usage{
find_targets(
  forward,
  reverse,
  assay_name,
  blast_path,
  blast_db,
  tax_db,
  threads,
  output_dir,
  banned_words,
  max_aligned_seqs,
  num_permutations,
  min_amp_len,
  max_amp_len
)
}
\arguments{
\item{forward}{forward primer}

\item{reverse}{reverse primer}

\item{assay_name}{name of the assay}

\item{blast_path}{path to blast version to use}

\item{blast_db}{path to blast database}

\item{tax_db}{path to taxonomy database created with make_tax_db()}

\item{threads}{number of threads to use}

\item{output_dir}{directory to write the files out to}

\item{banned_words}{list, a list of words used to identify sequences
with uncertain taxonomy}

\item{max_aligned_seqs}{Maximum number of sequences to keep for alignment
and data processing}

\item{num_permutations}{Maximum number of primer variants to include for
primers with ambiguous bases}

\item{min_amp_len}{minimum length of amplicon to be included in evaluation}

\item{max_amp_len}{maximum length of amplicon to be included in evaluation}
}
\value{
none
}
\description{
Run initial blast searches
}
