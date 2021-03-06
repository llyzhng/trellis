% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rearrangement-utils.R
\name{rearrangementType}
\alias{rearrangementType}
\title{Determine the type of rearrangement supported by each improper read pair}
\usage{
rearrangementType(object)
}
\arguments{
\item{object}{a \code{Rearrangement}}
}
\value{
a \code{data.frame} with colnames 'type' and 'percent'
}
\description{
Determine the type of rearrangement supported by each improper read pair
}
\details{
Rearrangements are typed by the following criteria for read 1 (R1)
  and read 2 (R2):

\enumerate{

  \item strand of R1 and R2
  \item orientation of R1 and R2 (e.g., R1 < R2)
  \item whether R1 and R2 are on different chromosomes
  \item whether R1 and R2 have an aberrant separation
}

\strong{Deletions:} A deletion results in a single new sequence
  junction.  There are two orientations of a R1 and R2 that support
  this junction: R1+ < R2- and R1- > R2+.  The former R1+ < R2- is
  the sequenced + strand DNA fragment, while the latter R1- > R2+
  is the sequenced - strand.  These orientations are the expected
  orientations in unrearranged genomes, except that the distance
  between R1 and its mate is larger than expected (> 10kb).

\strong{Amplicons:}
For purpose of downstream analyses, we only type intra-chromosomal
  amplicons that are replicated at the same site in the genome and
  without any changes to the strand orientation.  Illustration:

\preformatted{
 Reference:
        1  2  3     4 5 6    7 8 9
  +5' ------------|--------|-------------
  -3' ------------|--------|-------------

 Tumor
        1  2  3     4 5 6    4 5 6   7 8 9
  +5' ------------|--------X-------|---------
  -3' ------------|--------X-------|---------
 }
Note that only 'X' is a new sequence junction not seen in the
  reference.  Characteristics of 'X' are R1+ > R2- (+ fragment) and
  R1- < R2+ (- fragment).  The amplicon could also insert further
  downstream:

\preformatted{
 Tumor 
        1  2  3     4 5 6    7 8 9   4 5 6
  +5' ------------|--------|-------X---------
  -3' ------------|--------|-------X---------

 or further upstream:

 Tumor 
        4  5  6     1 2 3    4 5 6   7 8 9
  +5' ------------X--------|-------|---------
  -3' ------------X--------|-------|---------

}

  In each case, the rearrangement junctions given by 'X' would be
  identified by R1+ > R2 - and R1- < R2+.  

\strong{Inter-chromosomal translocations:} For translocations,
  positional orientation is determined by chromosome number and not
  genomic position.  WLOG, we consider chr1 to be less than chr2
  and chr22 to be less that chrX.  In the case of an unbalanced
  translocation, a single rearrangement junction will be identified
  by the improper read pairs.

\preformatted{
Reference   (-- chr1, == chr2)
  
     chr1  1 2 3   4 5 6   chr2  20 21 22   23  24  25
5' + ------------|------   ===============|==============
3' - ------------|------   ===============|==============

Tumor

      chr1  1 2 3   23 24 25 
5' + ------------|========== 
3' - ------------|==========

}

Read orientations R1+ < R2- and R1- > R2+ are consistent with the
  fusion of the two positive strands and the two negative strands,
  respectively.  If the translocation is balanced, we would also
  observe

\preformatted{

Tumor
     chr2  20 21 22   4 5 6
5' + ==============|--------
3' - ==============|--------

}

with R1+ > R2- and R1- < R2+.  Hence, there are four distinct read
  pair orientations for a balanced tranlocation. For the purpose of
  assessing gene fusions, we treat all translocations as if they
  are balanced even if we do not observe all 4 possible
  orientations.  An inversion translocation is typed and analyzed
  as an inversion.


\strong{Inversions:}

An intrachromosomal inversion: 

\preformatted{

Reference:  (-- positive strand,  == negative strand)

      1 2 3   4 5 6     7 8 9
5'+  ------>|-------->|------->
3'-  <======|<========|<=======

Tumor (inversion)

      1 2 3   6 5 4    7 8 9
5'+  ------>X=======>X------->
3'-  <======X<-------X<=======

}

Again, X's denote the new sequence junctions formed as a result of
  the inversion.  The left-most X's are supported by R1+ < R2+ (the
  top strand in the diagram) and R1+ > R1+ (bottom strand).  The
  right-most X's are supported by R1- < R2- (top strand) and R1- >
  R2- (bottom strand).

An inversion can also involve a translocation.
}
\note{
Type \code{amp1,amp3} would never be identified because the
  junction does not involve an aberrant separation between read
  pairs.
}
