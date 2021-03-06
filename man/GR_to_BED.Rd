% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{GR_to_BED}
\alias{GR_to_BED}
\title{GRanges to BED}
\usage{
GR_to_BED(GR, bed_fname = NULL, trackline = FALSE)
}
\arguments{
\item{GR}{GRanges-object.}

\item{bed_fname}{Path to where file will be saved. Must end in ".bed". Defaults to tempfile().}

\item{trackline}{logical: Whether a dummy trackline should be included in the .bed file}
}
\value{
Path to file.
}
\description{
Saves an GRanges-object as a BED-file and returns the file location.
}
\details{
Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
}
\seealso{
\code{\link{tempdir}} \code{\link{tempfile}}
}
\author{
Malte Thodberg
}
