% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusion-utils.R
\name{fusionTable}
\alias{fusionTable}
\title{Determine all possible fusions of a Rearrangement object}
\usage{
fusionTable(robj, txdb, tx, cds, genome, orgdb, id = "")
}
\arguments{
\item{robj}{A \code{Rearrangement} object}

\item{txdb}{A \code{TxDb} object}

\item{tx}{A \code{GRanges} object of transcripts from the \code{TxDb} object}

\item{cds}{A \code{GRangesList} object containing the CDS for each transcript}

\item{genome}{A \code{BSgenome} object}

\item{orgdb}{A \code{OrdDb} object for the Hsapiens genome}

\item{id}{A length-one character vector of the sample identifier}
}
\description{
Determine all possible fusions of a \code{Rearrangement} object and
evaluate whether each is in-frame.
}
\details{
Fusions are analyzed as follows.  First, all CDS from
  genes overlapping a rearrangement are extracted (ignoring strand)
  using the \code{getCDS} function.  In addition to extracting all
  the CDS, this function returns the possible gene fusion that may
  result from a rearrangement based on the modal rearrangement type
  (the modal rearrangement type is inferred from the strand and
  position orientation of read pairs).  The genes are denoted
  generically by

\preformatted{
            A          C
 5+ --------------|---------
 3- --------B-----|----D----

where "|" denotes a new sequence junction in a rearranged genome.
  } For a given fusion (say AC), we then \code{clip} CDS from A and
  CDS from B that are absent in the fused product.  After clipping,
  we \code{fuse} the remaining CDS from genes A and C. The function
  \code{tumorProtein} is used to derive the amino acid sequence of
  the tumor protein -- the protein that would be formed by fusing
  the tripped CDS from genes A and C.  To assess whether the fusion
  is in frame, we extract all known full transcripts from genes A
  and C and translate the DNA sequence of each transcript to an
  amino acid sequence. We refer to the amino acid sequences of the
  full CDS as the reference protein.  The function
  \code{referenceProtein} is a wrapper for getting the reference
  amino acid sequences.  Given the amino acid sequence of the
  clipped and fused transcripts (fused tumor protein) and the amino
  acid sequence of the full, unclipped transcripts (reference
  protein), we compare their sequences to assess whether the fusion
  is in-frame using the function \code{inFrameFusions}. The results
  are summarized in tabular format by the function
  \code{.fusionTable}.
}
\examples{
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.refGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
genome <- BSgenome.Hsapiens.UCSC.hg19
tx <- transcripts(txdb)
options(warn=-1)
cds.all <- cdsBy(txdb, "tx", use.names=TRUE)
data(rear_list)
r <- rear_list[["18557-18736"]]
fusionTable(r, txdb, tx, cds.all, genome,
            org.Hs.eg.db, id="test")
## in-frame fusion of CTNND2 and TRIO
}
\seealso{
See \code{\link{getCDS}} for how the CDS from genes
  involved in a rearrangement are extracted, \code{\link{clip}} and
  \code{\link{fuse}} for how transcripts are clipped and then
  fused, respectively.  See \code{\link{referenceProtein}} and
  \code{\link{tumorProtein}} for deriving germline (unrearranged)
  and somatic (rearranged) amino acid sequences.
}
