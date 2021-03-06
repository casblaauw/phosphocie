% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_kinases.R
\name{read_netphorest}
\alias{read_netphorest}
\title{Read in NetPhorest prediction output}
\usage{
read_netphorest(
  path,
  return_long = FALSE,
  cols_to_keep = c("fasta_id", "position", "residue", "fragment_11")
)
}
\arguments{
\item{path}{Path to the netphorest output file}

\item{return_long}{Optional boolean: whether to return the data in long format
(kinases on rows, like raw data), or whether to widen data (kinases as columns,
for use with U-CIE). Default is false (return wide)}

\item{cols_to_keep}{Columns to keep when pivoting to wide data.
By default, keeps 'fasta_id', 'position', 'residue', and 'fragment_11',
and drops 'method', 'organism', and 'binder_type'.}
}
\value{
A tibble with site/protein data and kinase scores
}
\description{
Read in NetPhorest input. Keeps only kinase scores.
}
\examples{
kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)

}
