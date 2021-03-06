% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepTaxdb.R
\name{prepare_tax_db}
\alias{prepare_tax_db}
\title{Prepare taxonomy database to use with eval_assay}
\usage{
prepare_tax_db(output_dir, tax_db)
}
\arguments{
\item{output_dir}{character, directory to output intermediate files}

\item{tax_db}{character, path to write taxonomy database}
}
\description{
This creates the taxonomy database used within bsPrimerTree.
This will download the data directly from NCBI and take roughly half an hour.
I recommend updating your blast database at the same time to avoid
inconsistencies between this database and the nt database.
}
\details{
To use this, you need the code available at
https://github.com/MVesuviusC/getTaxa
}
\examples{
\dontrun{
prepare_tax_db(output_dir = "\.", tax_db = "taxonomy.db")
}
}
