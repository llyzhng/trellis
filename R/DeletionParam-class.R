#' @include help.R
NULL

#' Parameter class for calling deletions
#'
#' @slot homozygous_thr length-one numeric vector
#'   indicating the maximum fold change expected for a homozygous
#'   deletion
#'
#' @slot hemizygous_thr length-one numeric vector
#'   indicating the maximum fold change expected for a hemizygous
#'   deletion
#'
#' @slot min_width length-one numeric vector indicating the minimum
#'   width of a CNV segment.
#'
#' @slot max_width length-one numeric vector indicating the maximum
#'   width of a CNV segment.
#'
#' @slot min_gapwidth length-one numeric vector used to reduce
#'   deletion calls. Provided that the segment means of adjacent
#'   deletions are qualitatively similar, two deletions less than
#'   min_gapwidth apart will be merged into a single CNV.
#'
#' @slot min_RPs_in_hemizygous length-one numeric vector indicating
#'   the minimum number of rearranged read pairs that can be spanned
#'   by a hemizygous deletion. [PARAMETER NOT USED]
#'
#' @slot max_RPs_in_homozygous length-one numeric vector indicating
#'   the maximum number of read pairs allowed to be spanned by a
#'   homozygous deletion.  [PARAMETER NOT USED]
#'
#' @slot nflanking_hemizygous length-one numeric vector indicating the
#'   minimum number of rearranged RPs required for calling a
#'   hemizygous deletion with higher confidence.  Hemizygous deletions
#'   with rearranged RPs are given higher priority over hemizygous
#'   deletions supported only by the tag density.
#'
#' @slot nflanking_homozygous length-one numeric vector indicating the
#'   minimum number of rearranged RPs required for calling a
#'   homozygous deletion with higher confidence.  Homozygous deletions
#'   with rearranged RPs are given higher priority over homozygous
#'   deletions supported only by the tag density.
#'
#' @slot maxRatioObservedToExpected length-one numeric vector [IGNORED?]
#'
#' @slot max_proportion_in_filter length-one numeric vector indicating
#'   the maximum proportion of a CNV that can be covered by one of the
#'   filters.  The filters are germline CNV, low mappability, low GC /
#'   unassembled regions, and outliers observed in a bin for 2 or more
#'   normal samples.
#'
#' @slot bam_seqlevels_style character string indicating style of seqnames
#'
#' @slot tumor_purity estimated fraction of sequenced DNA that is tumor-derived
#'
#' @slot remove_hemizygous whether to remove "hemizygous" calls ('hemizygous+' calls are not removed).
#'
#' @examples
#' DeletionParam()
#' @rdname DeletionParam-class
#' @export
setClass("DeletionParam", representation(homozygous_thr="numeric",
                                         hemizygous_thr="numeric",
                                         tumor_purity="numeric",
                                         min_width="numeric",
                                         max_width="numeric",
                                         min_gapwidth="numeric",
                                         min_RPs_in_hemizygous="numeric",
                                         max_RPs_in_homozygous="numeric",
                                         nflanking_hemizygous="integer",
                                         nflanking_homozygous="integer",
                                         maxRatioObservedToExpected="numeric",
                                         max_proportion_in_filter="numeric",
                                         bam_seqlevels_style="character",
                                         remove_hemizygous="logical"))


#' @param homozygous_thr length-one numeric vector indicating the maximum log2
#'   ratio expected for a homozygous deletion
#' @param hemizygous_thr length-one numeric vector indicating the maximum log2
#'   ratio expected for a hemizygous deletion
#' @param tumor_purity length-one numeric vector in [0, 1]
#' @param min_width length-one numeric vector indicating the minimum width of a
#'   CNV segment.
#' @param max_width length-one numeric vector indicating the maximum width of a
#'   CNV segment.
#' @param min_gapwidth length-one numeric vector used to reduce deletion calls.
#'   Provided that the segment means of adjacent deletions are qualitatively
#'   similar, two deletions less than min_gapwidth apart will be merged into a
#'   single CNV.
#' @param min_RPs_in_hemizygous length-one numeric vector indicating the minimum
#'   number of rearranged read pairs that can be spanned by a hemizygous
#'   deletion. [PARAMETER NOT USED]
#' @param max_RPs_in_homozygous length-one numeric vector indicating the maximum
#'   number of read pairs allowed to be spanned by a homozygous deletion.
#'   [PARAMETER NOT USED]
#' @param nflanking_hemizygous length-one numeric vector indicating the minimum
#'   number of rearranged RPs required for calling a hemizygous deletion with
#'   higher confidence. Hemizygous deletions with rearranged RPs are given
#'   higher priority over hemizygous deletions supported only by the tag
#'   density.
#' @param nflanking_homozygous length-one numeric vector indicating the minimum
#'   number of rearranged RPs required for calling a homozygous deletion with
#'   higher confidence. Homozygous deletions with rearranged RPs are given
#'   higher priority over homozygous deletions supported only by the tag
#'   density.
#' @param maxRatioObservedToExpected length-one numeric vector [IGNORED?]
#' @param max_proportion_in_filter length-one numeric vector indicating the
#'   maximum proportion of a CNV that can be covered by one of the filters. The
#'   filters are germline CNV, low mappability, low GC / unassembled regions,
#'   and outliers observed in a bin for 2 or more normal samples.
#' @param bam_seqlevels_style character string indicating style of seqnames
#' @param build length-1 character vector indicating genome build (must be hg18
#'   or hg19)
#' @export
#' @return a \code{DeletionParam} object
#' @rdname DeletionParam-class
DeletionParam <- function(min_width=2e3,
                          max_width=2e6,
                          min_gapwidth=10e3,
                          homozygous_thr=-3,
                          hemizygous_thr=-0.7,
                          tumor_purity=1,
                          min_RPs_in_hemizygous=1,
                          max_RPs_in_homozygous=100,
                          nflanking_hemizygous=5L,
                          nflanking_homozygous=5L,
                          maxRatioObservedToExpected=2L,
                          max_proportion_in_filter=0.75,
                          bam_seqlevels_style="UCSC",
                          remove_hemizygous=TRUE){
  new("DeletionParam",
      min_width=min_width,
      max_width=max_width,
      min_gapwidth=min_gapwidth,
      homozygous_thr=homozygous_thr,
      hemizygous_thr=hemizygous_thr,
      tumor_purity=tumor_purity, ## this should be estimated
      min_RPs_in_hemizygous=min_RPs_in_hemizygous,
      max_RPs_in_homozygous=max_RPs_in_homozygous,
      nflanking_hemizygous=nflanking_hemizygous,
      nflanking_homozygous=nflanking_homozygous,
      maxRatioObservedToExpected=maxRatioObservedToExpected,
      max_proportion_in_filter=max_proportion_in_filter,
      bam_seqlevels_style=bam_seqlevels_style,
      remove_hemizygous=remove_hemizygous)
}

#' @keywords internal
#' @rdname DeletionParam-class
#' @param object a \code{DeletionParam}
#' @export
bamSeqLevelsStyle <- function(object) object@bam_seqlevels_style

#' @keywords internal
#' @rdname DeletionParam-class
#' @param object a \code{DeletionParam}
#' @export
homozygousThr <- function(object) object@homozygous_thr

#' @keywords internal
#' @rdname DeletionParam-class
#' @param object a \code{DeletionParam}
#' @export
hemizygousThr <- function(object) object@hemizygous_thr

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
maximumWidth <- function(object) object@max_width

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
maximumProportionInFilter <- function(object) object@max_proportion_in_filter


#' @keywords internal
#' @rdname DeletionParam-class
#' @export
maxRatioObservedToExpected <- function(object) object@maxRatioObservedToExpected

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
minimumWidth <- function(object) object@min_width

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
minimumGapWidth <- function(object) object@min_gapwidth

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
minFlankingHemizygous <- function(object) object@nflanking_hemizygous

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
minFlankingHomozygous <- function(object) object@nflanking_homozygous

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
minRPsHemizygous <- function(object) object@min_RPs_in_hemizygous

#' @keywords internal
#' @rdname DeletionParam-class
#' @export
maxRPsHomozygous <- function(object) object@max_RPs_in_homozygous

setMethod("show", "DeletionParam", function(object){
  cat("'DeletionParam' class:\n")
  cat("   log_ratio threshold for homozygous: ", homozygousThr(object), "\n")
  cat("   log_ratio threshold for hemizygous: ", hemizygousThr(object), "\n")
  cat("   minimum # flanking read pairs\n")
  cat("            - homozygous:", minFlankingHomozygous(object), "\n")
  cat("            - hemizygous:", minFlankingHemizygous(object), "\n")
  cat("   maximum ratio of Observed/Expected:", maxRatioObservedToExpected(object), "\n")
})
