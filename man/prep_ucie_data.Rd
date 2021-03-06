% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ucie_transformations.R
\name{prep_ucie_data}
\alias{prep_ucie_data}
\title{Check and reformat data for use with U-CIE functions}
\usage{
prep_ucie_data(dataset, rownames_col = NULL, check_3D = TRUE)
}
\arguments{
\item{dataset}{A data frame or matrix of 2D or 3D values.}

\item{rownames_col}{Optional: name of column to ignore during calculations.}

\item{check_3D}{Optional: whether to check and force 3D output. Needed for reduced/cielab data,
not needed for KNN data.}
}
\description{
Takes a 2D or 3D matrix or data frame and returns a 3D matrix, checked for non-numeric values.
}
\details{
The function does a number of processing steps to accept as many data formats as feasible:
\itemize{
\item Rownames in the first column of the matrix (as tested by unconvertability to numeric) or of the data frame are discarded.
\item Character matrices are converted to numeric.
\item Data frames are converted to matrices.
\item If rownames_col is given, that column is discarded as well.
\item 2D matrices or data frames are expanded to 3D by padding with 1's.
}
}
\keyword{internal}
