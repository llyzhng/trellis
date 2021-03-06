% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rearlist-utils.R
\name{rFilters}
\alias{rFilters}
\title{Construct list of filters for somatic rearrangement analyses}
\usage{
rFilters(germline = GRanges(), rear = GRanges(),
  deletions = GRanges(), amplicons = GRanges())
}
\arguments{
\item{germline}{a \code{GRanges} object of germline filters}

\item{rear}{a \code{GRanges} object of rearrangements identified in germline samples (e.g., possibly artifactual)}

\item{deletions}{a \code{GRanges} object of deletions}

\item{amplicons}{a \code{GRanges} object of amplicons}
}
\value{
a \code{GRangesList} of rearrangement filter
}
\description{
Construct list of filters for somatic rearrangement analyses
}
