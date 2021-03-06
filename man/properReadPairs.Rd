% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment-utils.R
\name{properReadPairs}
\alias{properReadPairs}
\title{Import properly paired reads from a bam file}
\usage{
properReadPairs(bam_path, gr, param = DeletionParam())
}
\arguments{
\item{bam_path}{character string providing complete path to bam file}

\item{gr}{a \code{GRanges} object}

\item{param}{a \code{DeletionParam} object}
}
\value{
a \code{GAlignmentPairs} object
}
\description{
Import properly paired reads from a bam file
}
