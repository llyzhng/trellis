#' @include AllGenerics.R
NULL

#' Find total width of a \code{GRanges} object
#'
#' Adds the width of the intervals in a \code{GRanges} object. 
#'
#' @examples
#' gr <- GRanges(c("chr1", "chr2"), IRanges(c(11, 11), c(1000, 1000)))
#' totalWidth(gr)
#' @return numeric
#' @export
#' @param object a \code{GRanges} object
totalWidth <- function(object) sum(as.numeric(width(object)))


## Adapted from makeGAlignmentPairs in GenomicAlignments
makeGAlignmentPairs2 <- function(x, use.names=FALSE,
                                 use.mcols=FALSE, strandMode=1){
  if (!isTRUEorFALSE(use.names))
    stop("'use.names' must be TRUE or FALSE")
  if (!isTRUEorFALSE(use.mcols)) {
    if (!is.character(use.mcols))
      stop("'use.mcols' must be TRUE or FALSE or a character vector ",
           "specifying the metadata columns to propagate")
    if (!all(use.mcols %in% colnames(mcols(x))))
      stop("'use.mcols' must be a subset of 'colnames(mcols(x))'")
  }
  mate <- findMateAlignment(x)
  x_is_first <- GenomicAlignments:::.isFirstSegment.GAlignments(x)
  x_is_last <- GenomicAlignments:::.isLastSegment.GAlignments(x)
  first_idx <- which(!is.na(mate) & x_is_first)
  last_idx <- mate[first_idx]
  ## Fundamental property of the 'mate' vector: it's a permutation of order
  ## 2 and with no fixed point on the set of indices for which 'mate' is
  ## not NA.
  ## Check there are no fixed points.
  if (!all(first_idx != last_idx))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## Check order 2 (i.e. permuting a 2nd time brings back the original
  ## set of indices).
  if (!identical(mate[last_idx], first_idx))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## One more sanity check.
  if (!all(x_is_last[last_idx]))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## Check the 0x2 bit (isProperPair).
  x_is_proper <- as.logical(bamFlagAsBitMatrix(mcols(x)$flag,
                                               bitnames="isProperPair"))
  ans_is_proper <- x_is_proper[first_idx]
  ans_first <- x[first_idx]
  ans_last <- x[last_idx]
  ans_names <- NULL
  if (use.names)
    ans_names <- names(ans_first)
  names(ans_first) <- names(ans_last) <- NULL
  if (is.character(use.mcols)) {
    mcols(ans_first) <- mcols(ans_first)[use.mcols]
    mcols(ans_last) <- mcols(ans_last)[use.mcols]
  } else if (!use.mcols) {
    mcols(ans_first) <- mcols(ans_last) <- NULL
  }
  gps <- GAlignmentPairs(first=ans_first,
                         last=ans_last,
                         isProperPair=ans_is_proper,
                         names=ans_names,
                         strandMode=strandMode)
  gps
}

.trimInvalidReadsGAlign <- function(x){
  chromosome <- function(x) as.character(seqnames(x))
  ends <- setNames(end(x), chromosome(x))
  ends.info <- seqlengths(x)[names(ends)]
  is_valid <- ends <= ends.info
  x[is_valid]
}

#' Creates a default set of flags for reading improperly paired alignments
#'
#' These functions are wrappers for \code{scanBamFlag} and
#' \code{ScanBamParam} in the \code{Rsamtools} package.
#'
#' @seealso See \code{\link[Rsamtools]{scanBamFlag}} for complete details
#'   and \code{\link{improperAlignmentParams}} for a wrapper of this
#'   function that generates a \code{ScanBamParam} object using these
#'   flags.
#' 
#' @examples
#' require(Rsamtools)
#' 
#' flags <- scanBamFlag(isDuplicate=FALSE,
#'                      isPaired=TRUE,
#'                      isUnmappedQuery=FALSE,
#'                      hasUnmappedMate=FALSE,
#'                      isProperPair=FALSE)
#' flags2 <- improperAlignmentFlags()
#' identical(flags, flags2) #TRUE
#' print(flags)
#' @rdname alignment-flags
#' @export
improperAlignmentFlags <- function(){
  flags <- scanBamFlag(isDuplicate=FALSE,
                       isPaired=TRUE,
                       isUnmappedQuery=FALSE,
                       hasUnmappedMate=FALSE,
                       isProperPair=FALSE)
}


#' @rdname alignment-flags
#' @export
#' @examples
#'  
#' flags <- scanBamFlag(isDuplicate=FALSE,
#'                      isPaired=TRUE,
#'                      isUnmappedQuery=FALSE,
#'                      hasUnmappedMate=FALSE,
#'                      isProperPair=TRUE)
#' flags2 <- properAlignmentFlags()
#' identical(flags, flags2) #TRUE
#' print(flags)
properAlignmentFlags <- function(){
  flags <- scanBamFlag(isDuplicate=FALSE,
                       isPaired=TRUE,
                       isUnmappedQuery=FALSE,
                       hasUnmappedMate=FALSE,
                       isProperPair=TRUE)
}


#' @seealso See \code{\link[Rsamtools]{ScanBamParam}} and
#'   \code{\link[Rsamtools]{bamFlag}} in
#'   \code{Rsamtools} for full documentation.  See \code{improperAlignmentFlags} for the
#'   default set of flags.
#'
#' @examples
#'
#' flags <- improperAlignmentFlags()
#' print(flags)
#' params <- ScanBamParam(flag = flags, what=c("flag", "mrnm", "mpos", "mapq"))
#' params2 <- improperAlignmentParams()
#' print(params2)
#' identical(params, params2) #TRUE
#'
#' @param what A character vector (see \code{ScanBamParam} for details)
#' @param ... additional arguments to \code{ScanBamParam} such as \code{mapqFilter}
#' @param flag A length-two integer vector as provided by \code{improperAlignmentFlags}
#' @export
#' @rdname alignment-flags
improperAlignmentParams <- function(flag=improperAlignmentFlags(),
                                    what=c("flag", "mrnm", "mpos", "mapq",
                                           "qname"),
                                    ...){
  ScanBamParam(flag=flag, what=what, ...)
}


#' @export
#' @rdname alignment-flags
#' @examples 
#' flags <- properAlignmentFlags()
#' print(flags)
#' params <- ScanBamParam(flag = flags, what=c("flag", "mrnm", "mpos", "mapq"))
#' params2 <- properAlignmentParams()
#' identical(params, params2) #TRUE
#' print(params2)
properAlignmentParams <- function(flag=properAlignmentFlags(),
                                  what=c("flag", "mrnm", "mpos", "mapq"),
                                  ...){
  ScanBamParam(flag=flag, what=what, ...)
}


#' Extract all improperly paired reads from a bam file as a GAlignmentPairs object
#'
#' This function is a wrapper for readGAlignments. At the time this function was created, one could not create a GAlignmentPairs object of improperly paired reads using existing infrastructure in the GenomicAlignments package.  This function may be deprecated in the future.
#'
#' @return a \code{GAlignmentPairs} object
#' @keywords internal
#' @export
#' @param bam.file complete path to BAM file
#' @param param a \code{ScanBamParam} object.
#' @param build the reference genome buld that reads were aligned to.  Currently
#' supported builds include "hg19" and "hg18".
#'
#' @examples
#'   library(svbams)
#'   path <- system.file("extdata", package="svbams")
#'   bam.file <- file.path(path, "cgov10t.bam")
#'   irp <- getImproperAlignmentPairs(bam.file, build="hg19")
#'
#' @seealso See \code{\link[GenomicAlignments]{makeGAlignmentPairs}}
#'   for details regarding \code{use.mcols} argument.  See
#'   \code{\link{improperAlignmentParams}} for creating a
#'   \code{ScanBamParam} object with the appropriate flags for
#'   extracting improper read pairs.
getImproperAlignmentPairs <- function(bam.file,
                                      param=improperAlignmentParams(),
                                      build){
  flags <- improperAlignmentFlags()
  irp <- readGAlignments(bam.file, use.names=TRUE, param=param)
  irp <- .trimInvalidReadsGAlign(irp)
  irp2 <- makeGAlignmentPairs2(irp, use.mcols=TRUE, use.names=TRUE)
  genome(seqinfo(irp2)) <- build
  irp2
}

#' Extract properly paired reads from a bam file
#'
#' Function used by rearrangement analysis to extract the sequence of
#' properly paired reads from a bam file.
#'
#' @return a \code{GAlignmentPairs} object
#' @keywords internal
#' @export
#' @param bam.file a \code{BamViews} object
#' @param param a \code{ScanBamParam} object.
#' @param build the reference genome buld that reads were aligned to.  Currently
#' supported builds include "hg19" and "hg18".
#'
#' @examples 
#' library(svbams)
#' path <- system.file("extdata", package="svbams")
#' bam.file <- file.path(path, "cgov10t.bam")
#' irp <- getProperAlignmentPairs(bam.file, build="hg19")
#'
#' @seealso See \code{\link{properAlignmentParams}} for creating a
#'   \code{ScanBamParam} object with the appropriate flags for extracting
#'   proper read pairs.
getProperAlignmentPairs <- function(bam.file,
                                    param=properAlignmentParams(mapqFilter=0),
                                    build){
  ##bam.file <- bamPaths(object)
  irp <- readGAlignments(bam.file, use.names=TRUE, param=param)
  irp <- .trimInvalidReadsGAlign(irp)
  mapq_thr <- bamMapqFilter(param)
  irp2 <- makeGAlignmentPairs2(irp, use.mcols=TRUE, use.names=TRUE)
  genome(seqinfo(irp2)) <- build
  irp2
}

# Parse improper read pairs from a bam file
# 
# Parses improper read pairs from a bam file and saves the result as
# a serialized R object.  The file paths to the improper read pairs
# are given by accessors defined from the \code{AlignmentViews2}
# class.
# 
# @examples
#   library(Rsamtools)
#   require(TestBams)
#   extdata <- system.file("extdata", package="TestBams")
#   bam.file <- list.files(extdata, pattern="\\.bam$", full.name=TRUE)
#   bv <- BamViews(bam.file)
#   dp <- DataPaths(tempdir())
#   aviews <- AlignmentViews2(bv, dp)
#   \dontrun{
#     writeImproperAlignments2(aviews)
#     gps <- readRDS(file.path(dp["alignments/0improper"], rdsId(aviews)[1]))
#   }
# 
# @rdname AlignmentViews2
# @export
# @param bam.file complete path to BAM file
# @param param a \code{ScanBamParam} object
# @param mapq_thr  length-one numeric vector indicating lower limit of MAPQ score
# @param use.mcols length-one logical vector
# writeImproperAlignments2 <- function(bam.file,
#                                      param=improperAlignmentParams(mapqFilter=0)){
#   .Deprecated()
#   if(file.exists(bam.file)){
#     return(invisible())
#   }
#   ##aln_path <- improperPaths(aview)
#   irp <- getImproperAlignmentPairs(bam.file, param, build='hg19')
#   out.file <- improperPaths(aview)
#   saveRDS(irp, file=out.file)
#   invisible()
# }


thinReadPairQuery <- function(g, zoom.out=1){
  greater_20kb <- width(g) > 20e3
  g2 <- expandGRanges(g, zoom.out*width(g))
  ## for the big intervals, focus on the border areas
  if(!any(greater_20kb)){
    return(g2)
  }
  g3 <- g2[!greater_20kb]
  index <- greater_20kb
  g4 <- GRanges(seqnames(g)[index],
                IRanges(start(g)[greater_20kb]-5e3,
                        start(g)[greater_20kb]+5e3))
  g5 <- GRanges(seqnames(g)[index],
                IRanges(end(g)[greater_20kb]-5e3,
                        end(g)[greater_20kb]+5e3))
  gg <- c(granges(g3), g4, g5)
  dj <- disjoin(gg)
  ## TODO: consider reduce here. The same read pair will be in multiple bins
  dj
}


## will retrieve improper and proper reads
.mappedReadPairFlags <- function() scanBamFlag(isDuplicate=FALSE,
                                               isPaired=TRUE,
                                               isUnmappedQuery=FALSE,
                                               hasUnmappedMate=FALSE,
                                               isProperPair=TRUE)

.scan_all_readpairs <- function(granges, bam.file, flags){
  g <- reduce(granges)
  param <- ScanBamParam(flag=flags, what=c("flag", "mrnm", "mpos", "mapq"), which=g)
  ##x <- readGAlignmentsFromBam(bam.file, param=param, use.names=TRUE)
  x <- readGAlignments(bam.file, param=param, use.names=TRUE)
  x <- makeGAlignmentPairs2(x, use.mcols="flag")
  x
}

R1isFirst <- function(galp) start(first(galp)) < end(last(galp))

validFirstR1 <- function(galp){
  strand(first(galp)) == "+" & strand(last(galp)) == "-" & R1isFirst(galp)
}

validLastR1 <- function(galp){
  strand(first(galp)) == "-" & strand(last(galp)) == "+" & !R1isFirst(galp)
}

validPairForDeletion <- function(galp){
  ## for deletions:
  ## R1 +, R2-, R1<R2
  ## R1 -, R2+, R1>R2
  same_chrom <- chromosome(first(galp)) == chromosome(last(galp))
  ##if(!all(same_chrom)) galp <- galp[same_chrom]
  r1pos_r2neg <- validFirstR1(galp)
  r1neg_r2pos <- validLastR1(galp)
  r1pos_r2neg | r1neg_r2pos & same_chrom
}

#' Import properly paired reads from a bam file
#'
#'
#' @return a \code{GAlignmentPairs} object
#' @export
#' @param bam_path character string providing complete path to bam file
#' @param gr a \code{GRanges} object
#' @param param a \code{DeletionParam} object
properReadPairs <- function(bam_path, gr, param=DeletionParam()){
  query <- thinReadPairQuery(gr)
  if(length(query) == 0) {
    return(.GAlignmentPairs())
  }
  seqlevelsStyle(query) <- bamSeqLevelsStyle(param)
  flags <- .mappedReadPairFlags()
  galp <- .scan_all_readpairs(query,
                              bam.file=bam_path,
                              flags=flags)
  is_valid <- validPairForDeletion(galp)
  galp <- galp[is_valid]
  galp
}

#' Extract all mapped read pairs near a deletion(s)
#'
#' Wrapper for \code{readGAlignments} that extracts mapped read pairs
#' near deletions.
#'
#' @export
#' @return a \code{GAlignmentPairs} object
#' @param object a \code{GRanges} object
#' @param bam.file a character string providing complete path to bam file
readPairsNearVariant <- function(object, bam.file){
  flags <- .mappedReadPairFlags()
  ##granges <- thinReadPairQuery(object, thin)
  galp <- .scan_all_readpairs(object, bam.file=bam.file, flags=flags)
  seqlevelsStyle(galp) <- seqlevelsStyle(object)
  is_valid <- validPairForDeletion(galp)
  galp[is_valid]
}


#' Parse BAM file for improper read pairs near a set of GRanges
#'
#' All reads aligned to the intervals given by
#' \code{queryRanges(object)} are identified by the low-level function
#' \code{.scan_all_readpairs}.  This function reads alignments by
#' \code{readGAlignments} and then makes pairs of the alignments by
#' \code{makeGAlignmentPairs2}.  The latter function is an adaption of
#' the function \code{makeGAlignmentPairs} implemented in the
#' \code{GenomeAlignments} package but allows for the read pairs to be
#' improper.
#'
#' @param object Typically an \code{AmpliconGraph}, though the only
#'   requirement is that the method \code{queryRanges} is defined
#' @param bam.file character-vector providing valid complete path to a
#'   bam file
#' @param flags length-two integer vector as given by \code{scanBamFlags}
#' @return A \code{GAlignmentPairs} object 
#' 
#' @export
get_readpairs <- function(object, bam.file, flags=scanBamFlag()){
  g <- queryRanges(object)
  get_readpairs2(g, bam.file, flags)
}

#' Extract reads from a bam file
#' 
#' Parses a BAM file for read pairs near intervals specified by a \code{GRanges} object and 
#' returns the read pairs in a \code{GAlignmentPairs} object.
#' @param g A \code{GRanges} object
#' @param bam.file A character vector of the path to a bam file
#' @param flags A length-two integer vector as given by \code{scanBamFlags}
#' @seealso get_readpairs
#' @return A \code{GAlignmentPairs} object
#' @export
 get_readpairs2 <- function(g, bam.file, flags=scanBamFlag()){
   galp <- .scan_all_readpairs(g, bam.file, flags)
   validR1 <- overlapsAny(first(galp), g)
   validR2 <- overlapsAny(last(galp), g)
   galp <- galp[validR1 & validR2]
   galp
 }

#' Extract all improperly paired reads from an object with queryRanges defined
#'
#' This function is a wrapper for readGAlignments followed by
#' makeGAlignmentPairs2.  Only read pairs in which both the first read
#' and the last read in a pair overlap with the queryRanges are
#' returned.  Note, a given read pair need not overlap the same
#' queryRange.  To allow fuzzy matching of alignments to the
#' queryRanges, the queryRanges are expanded by 2kb in each direction.
#'
#' 
#' @export
#' 
#' @param object An object for which \code{queryRanges} method is defined
#' @param bam.file The complete file path to a bam file.
get_improper_readpairs <- function(object, bam.file){
  g <- reduce(queryRanges(object))
  if(totalWidth(g)==0) {
    galp <- GAlignmentPairs(first=GAlignments(), last=GAlignments(), isProperPair=logical())
    return(galp)
  }
  ## expand the query regions by 2kb on each side
  g2 <- reduce(expandGRanges(g, 2e3L))
  p <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE, isProperPair=FALSE),
                    what=c("flag", "mrnm", "mpos"), which=g)
  ##x <- readGAlignmentsFromBam(bam.file, param=p, use.names=TRUE)
  x <- readGAlignments(bam.file, param=p, use.names=TRUE)
  galp <- makeGAlignmentPairs2(x, use.mcols="flag")
  galp <- subsetByOverlaps2(galp, g)
  galp
}

subsetByOverlaps2 <- function(galp, g){
  validR1 <- overlapsAny(first(galp), g)
  validR2 <- overlapsAny(last(galp), g)
  galp <- galp[validR1 & validR2]
  galp
}

get_improper_readpairs2 <- function(g, bam.file){
  if(totalWidth(g)==0) {
    galp <- GAlignmentPairs(first=GAlignments(),
                            last=GAlignments(),
                            isProperPair=logical())
    return(galp)
  }
  ## expand the query regions by 2kb on each side
  g2 <- reduce(expandGRanges(g, 2e3L))
  p <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE, isProperPair=FALSE),
                    what=c("flag", "mrnm", "mpos"), which=g)
  ##x <- readGAlignmentsFromBam(bam.file, param=p, use.names=TRUE)
  x <- readGAlignments(bam.file, param=p, use.names=TRUE)
  galp <- makeGAlignmentPairs2(x, use.mcols="flag")
  validR1 <- overlapsAny(first(galp), g)
  validR2 <- overlapsAny(last(galp), g)
  galp <- galp[validR1 & validR2]
  galp
}

#' Convert GAlignmentPairs to GRanges while maintaining read pair information
#' 
#' Melts \code{GAlignmentPairs} objects to \code{GRanges} objects.
#'
#' @param ga a \code{GAlignments} object
#' @param is.improper a length-one logical vector.  \code{TRUE} if the reads are improperly paired, \code{FALSE} otherwise.
#' @param use.mcols logical for whether to keep the metadata columns of the \code{GAlignment} objects
#' 
#' @examples
#' library(svbams)
#' path <- system.file("extdata", package="svbams")
#' bam.file <- file.path(path, "cgov10t.bam")
#' irp <- getImproperAlignmentPairs(bam.file, build="hg19")
#' igr <- ga2gr(irp, is.improper=TRUE)
#' prp <- getProperAlignmentPairs(bam.file, build="hg19")
#' pgr <- ga2gr(prp, is.improper=FALSE)
#'
#' @return A \code{GRanges} object with metadata columns containing read pair information.
#'
#' @export
ga2gr <- function(ga, is.improper=FALSE, use.mcols=FALSE){
  id <- names(ga)
  r1.ga <- granges(first(ga), use.mcols=use.mcols)
  r2.ga <- granges(last(ga), use.mcols=use.mcols)
  names(r1.ga) <- names(r2.ga) <- NULL
  r1.ga$read <- rep("R1", length(r1.ga))
  r2.ga$read <- rep("R2", length(r2.ga))
  r1.ga$is.improper <- r2.ga$is.improper <- rep(is.improper, length(r1.ga))
  r1.ga$tagid <- id
  r2.ga$tagid <- id
  c(r1.ga, r2.ga)
}

#' Sort GRanges object with read pairs R1 and R2 by start of R1
#'
#' In order to plot the read pairs as position versus read pair index,
#' information on the mates needs to be kept intact and the visualization is
#' more clear if the reads are sorted by the first read in the pair. This
#' function turns the read tag in the 'id' field of the GRanges object to an
#' ordered factor. The levels of the factor are determined by the start position
#' of the first read (R1) in the pair.
#'
#' @param gr a \code{GRanges} object instantiated from a \code{GAlignmentPairs}
#' @return a \code{GRanges} object 
#' @examples
#'   library(svbams)
#'   library(TxDb.Hsapiens.UCSC.hg19.refGene)
#'   region <- GRanges("chr15", IRanges(63201003, 63209243))
#'   si <- seqinfo(TxDb.Hsapiens.UCSC.hg19.refGene)
#'   seqinfo(region) <- si["chr15", ]
#'
#'   path <-system.file("extdata", package="svbams")
#'   bampath <- list.files(path, pattern="cgov44t_revised.bam$",
#'                         full.names=TRUE)
#'
#'   iparams <- improperAlignmentParams(mapqFilter=30)
#'   pparams <- properAlignmentParams(mapqFilter=30)
#'   \dontrun{
#'     irp <- getImproperAlignmentPairs(bampath,
#'                                      iparams, build="hg19")
#'     g.irp <- ga2gr(irp, is.improper=TRUE)
#'     prp <- getProperAlignmentPairs(bampath,
#'                                    pparams, build="hg19")
#'     g.prp <- ga2gr(prp, is.improper=FALSE)
#'     gr <- c(g.irp, g.prp)
#'     gr <- sortByRead1(gr)
#'     gr
#' }
#' @export
sortByRead1 <- function(gr){
  gr.read1 <- gr[gr$read=="R1"]
  gr.read1 <- gr.read1[order(start(gr.read1))]
  tagid.levels <- as.character(gr.read1$tagid)
  gr$tagid <- factor(gr$tagid, levels=tagid.levels)
  gr <- gr[order(gr$tagid)]
  gr
}

#' Thin proper read pairs to reduce overplotting
#'
#' This function provides a simple interface to
#' subsample a \code{GRanges} object of properly paired
#' reads to reduce overplotting. 
#'
#' @param gr a \code{GRanges} object instantiated from a \code{GAlignmentPairs} object
#' @param thin integer indicating how much to thin the properly paired reads.
#' @details Setting the parameter \code{thin} to 10 (default) will 
#' return a \code{GRanges} object with 1/10 the original number of 
#' properly paired reads in \code{gr}.
#' @return A \code{GRanges} object  
#' @examples 
#' library(svbams)
#' path <- system.file("extdata", package="svbams")
#' bam.file <- file.path(path, "cgov10t.bam")
#' prp <- getProperAlignmentPairs(bam.file, build="hg19")
#' pgr <- ga2gr(prp, is.improper=FALSE)
#' length(pgr)
#' pgr2 <- thinProperPairs(pgr, 100)
#' length(pgr2)
#' @export
thinProperPairs <- function(gr, thin=10){
  is.proper <- !gr$is.improper
  if(!any(is.proper)){
    return(gr)
  }
  gr.prp <- gr[is.proper]
  gr.prp <- sortByRead1(gr.prp)
  tagid <- unique(as.character(gr.prp$tagid))
  tagid <- tagid[seq(1, length(tagid), thin)]
  gr.prp <- gr.prp[gr.prp$tagid %in% tagid]
  gr.prp
}

.seqdataframe <- function(gpairs, MAX){
  df <- as.data.frame(first(gpairs))
  if(nrow(df) > MAX){
    index <- sample(seq_len(MAX), MAX)
  } else index <- seq_len(nrow(df))
  df <- df[index, ]
  df2 <- as.data.frame(last(gpairs))
  df2 <- df2[index, ]
  df$read <- rep("R1", nrow(df))
  df2$read <- rep("R2", nrow(df))
  df <- rbind(df, df2)
  df
}

#' Assesses whether two reads in a pair are aberrantly separated with respect to the reference genome
#'
#' Determines whether separation between first and last read in a
#' GAlignmentPairs is greater than some distance.
#'
#' Evaluates to TRUE if
#' @param gpairs a \code{GAlignmentPairs} object
#' @param distance the minimum distance in base pairs between two reads in 
#' order to call them aberrantly separated
#' @return logical vector of the same length as \code{gpairs}
#' @export
aberrantSep <- function(gpairs, distance=10e3){
  abs((start(first(gpairs)) - start(last(gpairs)))) > distance
}

#' Assesses whether any read pairs are duplicates in a GAlignmentPairs object
#'
#' @param aln.pairs a \code{GAlignmentPairs} object
#' @return logical vector of same length as the \code{aln.pairs}
#' @export
isDuplicate <- function(aln.pairs){
  chr1 <- chromosome(first(aln.pairs))
  chr2 <- chromosome(last(aln.pairs))
  x1 <- start(first(aln.pairs))
  x2 <- end(last(aln.pairs))
  string <- paste(chr1, x1, chr2, x2, sep="_")
  duplicated(string)
}

#' Remove duplicates paired r
#' @export
#' @param gpairs TODO
#' @param bins TODO
#' @param params TODO
filterPairedReads <- function(gpairs, bins, params){
  gpairs <- gpairs[ aberrantSep(gpairs, rpSeparation(params)) ]
  gpairs <- gpairs [ !isDuplicate(gpairs) ]
  cnt <- countOverlaps(bins, first(gpairs)) + countOverlaps(bins, last(gpairs))
  gpairs <- subsetByOverlaps(gpairs, bins [ cnt >= minNumberTagsPerCluster(params) ])
  gpairs
}

.getReadSeqsForRear <- function(object, bam.file, param, MAX, build){
  irp <- improper(object)
  flags <- improperAlignmentFlags()
  lb <- linkedBins(object)
  bins <- c(granges(lb), granges(lb$linked.to))
  what <- c("qname", "rname", "seq", "flag", "mrnm", "mpos", "mapq")
  scan.param <- ScanBamParam(flag=flags, what=what, which=bins, mapqFilter=30)
  rps <- getImproperAlignmentPairs(bam.file,
                                   scan.param, 
                                   build = build)
  rps2 <- filterPairedReads(rps, bins, param)
  df <- .seqdataframe(rps2, MAX)
  ##df$id <- rep(names(align_view), nrow(df))
  df$id <- rep(basename(bam.file), nrow(df))
  df$rearrangement.id <- rep(names(linkedBins(object)), nrow(df))
  df
}

##.getReadSeqsForRear2 <- function(irp, bam.file, param, MAX){
##  flags <- improperAlignmentFlags()
##  ##lb <- linkedBins(object)
##  bins <- c(granges(lb), granges(lb$linked.to))
##  what <- c("qname", "rname", "seq", "flag", "mrnm", "mpos", "mapq")
##  scan.param <- ScanBamParam(flag=flags, what=what, which=bins, mapqFilter=30)
##  rps <- getImproperAlignmentPairs(bam.file,
##                                   scan.param)
##  rps2 <- filterPairedReads(rps, bins, param)
##  df <- .seqdataframe(rps2, MAX)
##  ##df$id <- rep(names(align_view), nrow(df))
##  df$id <- rep(basename(bam.file), nrow(df))
##  df$rearrangement.id <- rep(names(linkedBins(object)), nrow(df))
##  df
##}


#' @aliases getSequenceOfReads,Rearrangement-method
#' @rdname getSequenceOfReads-methods
setMethod("getSequenceOfReads", "Rearrangement",
          function(object,
                   bam.file,
                   params=RearrangementParams(),
                   MAX=25L, 
                   build){
            dat <- .getReadSeqsForRear(object, bam.file, params, MAX, build = build)
            dat2 <- as.tibble(dat)
            dat2
          })

#' @aliases getSequenceOfReads,RearrangementList-method
#' @rdname getSequenceOfReads-methods
setMethod("getSequenceOfReads", "RearrangementList",
          function(object,
                   bam.file,
                   params=RearrangementParams(),
                   MAX=25L, 
                   build){
            tag.list <- vector("list", length(object))
            for(j in seq_along(object)){
              r <- object[[j]]
              tag.list[[j]] <- .getReadSeqsForRear(r, bam.file,
                                                   params, MAX=MAX, 
                                                   build = build)
            }
            tags <- do.call(rbind, tag.list)
            readId <- paste(tags$qname, tags$read, sep="_")
            tags <- tags[!duplicated(readId), ]
            tags2 <- as.tibble(tags)
            tags2
          })
