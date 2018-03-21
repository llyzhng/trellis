#' Read output from the command-line blat
#'
#' The command-line version of blat can be downloaded from sourceforge:
#' \url{http://sourceforge.net/projects/blat/files/}
#'
#' @export
#' @param filename character string providing full path to blat output
#' @return a \code{data.frame} of alignment records from blat
#'
#' @seealso See \code{\link{blatScores}} for evaluating whether the
#'   blat alignments support a novel sequence junction
#'
readBlat <- function(filename){
  blat <- read.delim(filename, skip=2, nrows=5, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  nms <- paste0(as.character(blat[1, ]), as.character(blat[2, ]))
  nms <- gsub(" ", "", nms)
  ##nms[22:23] <- c("seq1", "seq2")
  blat <- read.delim(filename, skip=5, stringsAsFactors=FALSE, header=FALSE, sep="\t")
  colnames(blat) <- nms
  blat$Tname <- gsub(".fa", "", blat$Tname)
  blat
}

blatGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
               strand=factor(blat$strand, levels=c("+", "-", "*")))
  seqlevels(g, pruning.mode="coarse") <- sl
  g
}

elandGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$eland.chr, IRanges(blat$eland.start, blat$eland.end))
  seqlevels(g, pruning.mode="coarse") <- sl
  g
}

elandAndBlatOverlap <- function(blat){
  blat$Qname <- factor(blat$Qname, levels=unique(blat$Qname))
  g.blat <- blatGRanges(blat)
  g.eland <- elandGRanges(blat)
  if(length(g.blat) != length(g.eland)){
    msg <- "Check whether unique(blat$Tname) is the same as unique(blat$eland.chr)"
    stop(msg)
  }
  strand(g.eland) <- strand(g.blat)
  same_seqname <- chromosome(g.eland) == chromosome(g.blat)
  tmp <- pintersect(g.eland[same_seqname], g.blat[same_seqname]) ##resolve.empty="start.x")
  is_overlap <- rep(FALSE, length(g.eland))
  is_overlap[same_seqname] <- width(tmp) > 0
  is_overlap
}

annotateBlatRecords <- function(blat, tag.sequences){
  tmp <- tag.sequences
  qname2 <- paste0(tmp$qname, "_", tmp$read)
  query.sequences <- setNames(tmp$seq, qname2)
  query.sequences <- query.sequences[names(query.sequences) %in% blat$Qname]
  query.sequences <- query.sequences[blat$Qname]
  ##identical(names(query.sequences), blat$Qname)
  blat$Qsequence <- query.sequences

  qnames <- factor(blat$Qname, levels=unique(blat$Qname))
  indices <- split(1:nrow(blat), qnames)
  ## T is for target (reference genome)
  ## Q is for query (tag sequence)
  sampleids <- tag.sequences$id
  ##sampleids <- rep(sampleids, elementNROWS(tag.sequences))
  eland.starts <- tag.sequences$start
  eland.ends <- tag.sequences$end
  eland.chr <- tag.sequences$seqnames
  eland.strand <- tag.sequences$strand
  ##
  ## This no longer applies
  ##
  tagnames <- paste(tag.sequences$qname, tag.sequences$read, sep="_")
  names(sampleids) <- tagnames
  names(eland.starts) <- tagnames
  names(eland.ends) <- tagnames
  names(eland.chr) <- tagnames
  names(eland.strand) <- tagnames
  ##
  ##
  ##
  ##table(blat$Qname %in% tagnames)
  qname_in_tgs <- blat$Qname %in% tagnames
  if(!all(qname_in_tgs)){
    stop("Qnames in blat records not in tag names. One of the objects must be outdated.")
    blat <- blat[blat$Qname %in% tagnames, ]
  }
  blat$id <- sampleids[blat$Qname]
  blat$eland.chr <- eland.chr[blat$Qname]
  blat$eland.start <- eland.starts[blat$Qname]
  blat$eland.end <- eland.ends[blat$Qname]
  blat$eland.strand <- eland.strand[blat$Qname]
  ## replace Tsize with T.end-T.start
  blat$Tsize <- blat$Tend-blat$Tstart
  blat$is_overlap <- elandAndBlatOverlap(blat)
  rownames(blat) <- NULL
  blat
}

## A read can not have 2 hits meeting the 90-90 rule
removeReadsWithoutMate <- function(blat){
  rid <- blat$Qname
  rpid <- gsub("_R[12]$", "", rid)
  rpid <- factor(rpid, levels=unique(rpid))
  rplist <- lapply(split(rid, rpid), function(x) unique(x))
  el <- elementNROWS(rplist)
  ## there are 36 reads with no mate listed
  rps.without.mate <- names(rplist)[el != 2]
  blat <- blat[!rpid %in% rps.without.mate, ]
}

.listTagsByGroup <- function(tags, tag.group){
  tag.group <- factor(tag.group, levels=unique(tag.group))
  tagnames <- paste(tags$qname, tags$read, sep="_")
  tag.list <- split(tagnames, tag.group)
}

addXCoordinateForTag <- function(blat){
  qnm <- blat$Qname
  rid <- gsub("_R[12]", "", qnm)
  x <- as.integer(factor(rid, levels=unique(rid)))
  is_r2 <- rep(FALSE, length(x))
  is_r2[grep("_R2$", qnm)] <- TRUE
  ## for plotting
  x[!is_r2] <- x[!is_r2] - 0.2
  x[is_r2] <- x[is_r2] + 0.2
  x
}

blatStatsPerTag <- function(blat.records, tag_length){
  ## summary statistics for a single read
  is_90 <- blat.records$match > 90
  ## Tsize in blat.records is size of target sequence (chromosome)
  Tsize <- abs(blat.records$Tend-blat.records$Tstart)
  is_size_near100 <- Tsize > (tag_length - 1/5*tag_length) &
    Tsize < (tag_length + 1/5*tag_length)
  is_90 <- is_90 & is_size_near100
  is_overlap <- blat.records$is_overlap
  ## calculate number of near-perfect matches for each tag
  n.matches <- sapply(split(is_90, blat.records$Qname), sum)
  n.eland.matches <- sapply(split(is_overlap & is_90, blat.records$Qname), sum)
  ## there should only be one eland alignment with high quality (n.matches == 1)
  ## AND this read should overlap the eland alignment
  n.matches==1 & n.eland.matches==1
}

.blatStatsRearrangement <- function(blat, thr=0.8, tag_length){
  cols <- c("Qname", "match", "is_overlap", "Tstart", "Tend")
  blat <- blat[, cols]
  if(nrow(blat) == 0) return(NULL)
  blat$tag_index <- addXCoordinateForTag(blat)
  stats <- blatStatsPerTag(blat, tag_length)
  proportion_pass <- mean(stats)
  qnames <- gsub("R[12]", "", blat$Qname)
  n_tags <- length(unique(qnames))
  is_pass <- proportion_pass > thr & n_tags >= 5
  ##blat$passQC <- rep(proportion_pass > thr, nrow(blat))
  blat$passQC <- is_pass
  blat
}

#' blatScores assesses whether the improper read pairs at a
#' rearrangement junction provide strong support of the rearrangement
#'
#' In the following, we refer to a 'record' as one row in the table of
#' blat output -- i.e., one (of possibly many) alignments for a read.
#' This function adds an indicator for whether the blat alignment is
#' consistent with the original alignment (\code{is_overlap}), an
#' indicator of whether it passes quality control (\code{passQC}) (see
#' details), the id of the rearrangement (\code{rearrangement}), and
#' the sample id (\code{id}).
#'
#' @details
#'
#' \strong{Read-level QC:} 
#'
#' For each record, we evaluate
#'
#' \enumerate{
#' 
#' \item whether the matching score is at least 90
#' 
#' \item whether the size of the target alignment (Tstart - Tend) is
#'   less than 120bp and more than 80bp
#'
#' \item whether the blat alignment overlaps with any of the original
#'   alignment
#'
#' }
#'
#' The records are then grouped by Qname.  For each Qname, we compute
#'   the number of reads with a score above 90 and within the
#'   specified size range.  In addition, we compute the number of
#'   reads with a score above 90, within the specified size range, and
#'   that overlap with the original alignment.  The Qname passes QC if
#'   both sums evaluate to 1.  That is, a read passes QC only if there
#'   is a single record with a high BLAT score within the specified
#'   size range and this single record overlaps with the original
#'   alignment.
#'
#' \strong{Rearrangement-level QC:}
#'
#' Given a pass / fail designation for each read by the above
#'   analysis, we group the reads by the rearrangement id.  A
#'   rearrangement passes QC if > 80% of the reads supporting the
#'   rearrangement pass QC.
#' 
#' @examples
#'  qnames <- c(paste0(letters[1:10], "_R1"),
#'              paste0(letters[1:10], "_R2"))
#'  ## only 1 location
#'  numberAlignedLocations <- rep(1, length(qnames))
#'  matchScores <- rep(95, length(qnames))
#'  Tend <-  150
#'  Tstart <- 51
#'  ## Made up output from blat
#'  sblat <- data.frame(Qname=qnames,
#'                      match=matchScores,
#'                      Tend=Tend,
#'                      Tstart=Tstart,
#'                      Tname=rep("chr1", length(qnames)),
#'                      strand=rep("+", length(qnames)),
#'                      stringsAsFactors=FALSE)
#' 
#' ## Made up information on a rearrangement, including the locations
#' ##  at which the reads were originally aligned
#' 
#'  stags <- data.frame(qname=rep(letters[1:10], 2),
#'                      read=rep(c("R1", "R2"), each=10),
#'                      seqnames=rep("chr1", 10),
#'                      strand=rep("+", 10),
#'                      start=Tstart,
#'                      end=Tend,
#'                      seq=replicate(10, paste(sample(c("g","c"), 10, replace=TRUE),
#'                          collapse="")),
#'                      stringsAsFactors=FALSE)
#'  stags$id <- "CGOV32T"
#'  stags$rearrangement.id <- "1-3"
#'  ## For each tag, calculate the number of near-perfect matches of the
#'  ## right size that overlap with eland.  If the number is 0 or more
#'  ## than 1, then the tag 'fails'.
#'  blat <- svalignments:::annotateBlatRecords(sblat, stags)
#'  nrow(blat) == 20
#'  s <- blatScores(sblat, stags, "SOME_ID")
#'  head(s)
#' @export
#' @param blat a \code{data.frame} of results from  command-line blat
#' @param tags a \code{data.frame} containing read names and the original alignment locations
#' @param id  a character-vector of sample identifiers
#' @param thr a length-one numeric vector indicating the fraction of
#'   reads at a rearrangement that must pass the read-level QC.
blatScores <- function(blat, tags, id, thr=0.8){
  tag_length <- nchar(tags$seq[1])
  blat$match <- blat$match/tag_length * 100
  rownames(tags) <- paste0(tags$qname, "_", tags$read)
  blat <- annotateBlatRecords(blat, tags)
  blat <- removeReadsWithoutMate(blat)
  tagnames.list <- .listTagsByGroup(tags, tags[["rearrangement.id"]])
  blat.parsed <- vector("list", length(tagnames.list))
  for(j in seq_along(tagnames.list)){
    tagnames <- tagnames.list[[j]]
    rid <- names(tagnames.list)[j]
    blat_rid <- blat[blat$Qname %in% tagnames, ]
    result <- .blatStatsRearrangement(blat_rid, thr=thr, tag_length)
    if(!is.null(result)){
      result$rearrangement <- rep(rid, nrow(result))
    }
    blat.parsed[[j]] <- result
  }
  blat.parsed <- blat.parsed[!sapply(blat.parsed, is.null)]
  blat.parsed <- do.call("rbind", blat.parsed)
  blat.parsed$id <- rep(id, nrow(blat.parsed))
  blat.parsed
}


#' Reads files output by the blat executable.
#'
#' The blat executable is run with the default set of arguments using
#' a reference genome comprised of standard sequences without
#' alternates.  This function is a wrapper for \code{blatScores}. If a
#' file containing \code{blatScores} already exists, this function
#' only reads previously saved computations from disk.
#'
#' @seealso See \code{blatScores}
#' @examples
#' \dontrun{
#'  library(svovarian)
#'  dirs <- projectOvarian(rootname="OvarianData2")
#'  if(FALSE){
#'    tags <- readRDS(file.path(dirs[["3read"]], "CGOV2T.rds"))
#'    blat <- readBlat(file.path(dirs["1blat"], "CGOV2T.txt"))
#'    parsed <- scoreBlatExperiment(id, blat, tags, dirs)
#'    parsed
#'  }
#' }
#'
#' @return a list of \code{tbl_df} objects
#'
#' @param id sample id
#'
#' @param blat data.frame of blat records from the command-line blat tool
#'
#' @param tags a \code{tbl_df} object of the read sequences
#'
#' @param dirs a \code{DataPaths} object
#'
#' @param thr a length-one numeric vector indicating the fraction of
#'   improper reads that are of high quality for a given
#'   rearrangement.
#'
scoreBlatExperiment <- function(id, blat, tags, dirs, thr=0.8){
  parsed_file <- file.path(dirs[["4parsed_mapped"]], paste0(id, ".rds"))
  if(file.exists(parsed_file)){
    blat2 <- readRDS(parsed_file)
  } else {
    blat2 <- blatScores(blat, tags, id=id, thr=thr)
    saveRDS(blat2, file=parsed_file)
  }
  blat2 <- tbl_df(blat2)
  blat2
}

overlapsBlatRecord <- function(linked_bins, blat_record, maxgap=200){
  overlapsAny(linked_bins, blat_record, maxgap=maxgap) &
    overlapsAny(linked_bins$linked.to, blat_record, maxgap=maxgap)
}

sequenceRanges <- function(blat){
  ##
  ## We ignore the seqname, but we use GRanges instead of IRanges to
  ## make use of the metadata (here, the score from blat)
  ##
  GRanges("seq", IRanges(blat$qstart, blat$qend), match=blat$match)
}

sequenceRanges2 <- function(blat){
  ##
  ## We ignore the seqname, but we use GRanges instead of IRanges to
  ## make use of the metadata (here, the score from blat)
  ##
  GRanges("seq", IRanges(blat$qStarts, width=blat$blockSizes),
          match=blat$match,
          Qsize=blat$Qsize)
}

multipleAlignmentRecords <- function(records){
  records <- records[ records$match < 95]
  recordlist <- split(records, records$qname)
  n.alignments <- elementNROWS(recordlist)
  recordlist <- recordlist[n.alignments >= 2]
  unlist(GRangesList(recordlist))
}

integer_vector <- function(x){
  as.integer(unlist(strsplit(x, ",")))
}

##start_vector <- function(g$tStarts){
##  as.integer(unlist(strsplit(tStarts, ",")))
##}
##
##block_vector <- function(blockSizes){
##  as.integer(unlist(strsplit(blockSizes, ",")))
##}

blatStartList <- function(blat.gr){
  lapply(blat.gr$tStarts, integer_vector)
}

blatBlockList <- function(blat.gr){
  lapply(blat.gr$blocksizes, integer_vector)
}

candidateSplitRead <- function(blat.gr){
  blat.gr$match < 95 | blat.gr$blockcount > 1
}

numberAlignmentRecords <- function(blat.gr){
  L <- length(blat.gr)
  if(L == 1){
    L <- blat.gr$blockcount
  }
  L
}

.each_block_granges <- function(g){
  starts <- integer_vector(g$tStarts)
  L <- length(starts)
  widths <- integer_vector(g$blockSizes)
  qstarts <- integer_vector(g$qStarts)
  bsizes <- integer_vector(g$blockSizes)
  chrom <- rep(chromosome(g), g$blockcount)
  bmatch <- rep(g$match, g$blockcount)
  gapbases <- rep(g$gapbases, g$blockcount)
  g2 <- GRanges(chrom,
                IRanges(starts, width=widths),
                qStarts=qstarts,
                blockSizes=bsizes,
                gapbases=gapbases,
                match=bmatch)
  gr <- reduce(g2, with.revmap=TRUE)
  revmap <- mcols(gr)$revmap
  tmp <- relist(g2$qStarts[unlist(revmap)], revmap)
  qstarts <- sapply(tmp, min)
  tmp <- relist(g2$blockSizes[unlist(revmap)], revmap)
  bsizes <- sapply(tmp, sum)
  tmp <- relist(g2$match[unlist(revmap)], revmap)
  bmatch <- sapply(tmp, mean)
  tmp <- relist(g2$gapbases[unlist(revmap)], revmap)
  gapbases <- sapply(tmp, min)
  L <- length(gr)
  gr$qStarts <- qstarts
  gr$blockSizes <- bsizes
  gr$match <- bmatch
  gr$rear.id <- rep(g$rear.id[1], L)
  gr$qname <- rep(g$qname[1], L)
  gr$gapbases <- gapbases
  gr$Qsize <- rep(g$Qsize[1], L)
  gr
}

eachBlockAsGRanges <- function(blat.grl){
##  browser()
##  for(i in seq_along(blat.grl)){
##    .each_block_granges(blat.grl[[i]])
##  }
  blat.grl2 <- lapply(blat.grl, .each_block_granges)
  GRangesList(blat.grl2)
}

multipleAlignmentRecords2 <- function(records){
  records <- records[ candidateSplitRead(records) ]
  recordlist <- split(records, records$qname)
  ##n.alignments <- elementNROWS(recordlist)
  n.alignments <- sapply(recordlist, numberAlignmentRecords)
  recordlist <- recordlist[n.alignments >= 2]
  unlist(GRangesList(recordlist))
}

.splitread_intersection <- function(x){
  w <- width(intersect(x[1], x[2]))
  ## return 0 if no overlap
  if(length(w) == 0) w <- 0L 
  w
}

splitreadIntersection <- function(g){
  sapply(g, .splitread_intersection)
}

#' Identify rearranged reads -- initiallly unmapped reads that can be
#' aligned by blat to span a novel sequence junction.
#'
#' @return a `GRangesList` of blat records that map to both sides of a
#' sequence  junction. Each list element corresponds to one read that is
#'  aligned to two  locations (i.e., each element of the list consists of the
#'  vector of reads that supports one rearrangement).
#'
#'
#'
#' @export
#'
#' @param linked_bins a \code{GRanges} of linked bins (e.g., gotten by \code{linkedBins(rearrangement.list)})
#'
#' @param blat a data.frame of blat alignment records
#'
#' @param maxgap this maximum gap between the mapped read and the
#'   genomic intervals of the improper read clusters
rearrangedReads <- function(linked_bins, blat, maxgap=500){
  ##lb <- linkedBins(rlist)
  ## filter blat alignments
  is_na <- is.na(blat$Tstart)
  if(any(is_na)){
    blat <- blat[!is_na, ]
  }
  blat <- blat[blat$Tname %in% seqlevels(linked_bins), ]
  blat_gr <- blat_to_granges(blat, linked_bins)
  genome(blat_gr) <- genome(linked_bins)
  ##
  ## A blat record must overlap one of the intervals in a linked bin
  ##
  is.overlap <- overlapsLinkedBins(blat_gr, linked_bins, maxgap=maxgap)
  blat_gr <- blat_gr[ is.overlap ]
  ## We are looking for split read alignments--a read must have at
  ## least 2 alignments.  Get rid of all reads with only a single
  ## alignment, and all alignments that have a match score greater
  ## than 95.
  blat_gr <- multipleAlignmentRecords2(blat_gr)
  ##
  ## Both intervals in a linked bin must overlap a blat record
  ##
  overlaps_both <- overlapsBlatRecord(linked_bins, blat_gr, maxgap)
  sr_list <- vector("list", length(linked_bins))
  names(sr_list) <- names(linked_bins)
  for(i in seq_along(linked_bins)){
    if(!overlaps_both[i]) sr_list[[i]] <- empty_record()
    sr_list[[i]] <-   rearrangedReads2(linked_bins[i], blat_gr, maxgap)
  }
  sr.grl <- GRangesList(sr_list)
  sr.grl
}

empty_record <- function(){
  g <- GRanges()
  g$revmap <- IntegerList()
  g$qStarts <- integer()
  g$blockSizes <- integer()
  g$match <- numeric()
  g$rear.id <- character()
  g$qname <- character()
  g$gapbases <- integer()
  g$Qsize <- integer()
  g
}

blat_to_granges <- function(blat, lb){
  GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
          match=blat$match, qname=blat$Qname,
          qstart=blat$Qstart,
          qend=blat$Qend,
          tStarts=blat$tStarts,
          blockSizes=blat$blockSizes,
          gapbases=blat$Tgapbases,
          blockcount=blat$blockcount,
          qStarts=blat$qStarts,
          Qsize=blat$Qsize,
          seqinfo=seqinfo(lb))
}

overlapsLinkedBins <- function(blat_gr, lb, maxgap=500){
  overlapsAny(blat_gr, lb, maxgap=maxgap) |
    overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)  
}

get_rearrangement_id <- function(records, lb, maxgap=500){
  records$rear.id <- NA
  h1 <- findOverlaps(records, lb, maxgap=maxgap)
  records$rear.id[queryHits(h1)] <- names(lb)[subjectHits(h1)]
  h2 <- findOverlaps(records, lb$linked.to, maxgap=maxgap)
  records$rear.id[queryHits(h2)] <- names(lb)[subjectHits(h2)]
  records$rear.id
}


rearrangedReads2 <- function(linked_bins, blat_gr, maxgap=500){
  lb <- linked_bins
  if(is.null(names(lb))) stop("Expect linked bins to be named by rearrangement id")
  record_overlaps <- overlapsAny(blat_gr, lb, maxgap=maxgap) |
    overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)
  records <- blat_gr[ record_overlaps ]
  if(length(records) == 0){
    return(GRanges())
  }
  ##
  ## Assign each sequence (qname) to a rearrangement
  ##
  records$rear.id <- get_rearrangement_id(records, lb)
  ##
  ## split records by qname
  ##
  records_qname <- split(records, records$qname)
  ##
  ## Each sequence should have a unique rearrangement
  ##
  n.rear <- sapply(records_qname, function(x) length(unique(x$rear.id)))
  records_qname <- records_qname [ n.rear == 1 ]
  ##
  ## And the number of records for a given rearrangement should be 2
  ##
  blat.grl <- eachBlockAsGRanges(records_qname)

  # Pass along seqinfo
  seqlevels(blat.grl) <- seqlevels(lb)
  seqlengths(blat.grl) <- seqlengths(lb)
  genome(blat.grl) <- genome(lb)

  overlapFun <- function(g1, lb, ...){
    o1 <- any(overlapsAny(g1, lb, ...))
    o2 <- any(overlapsAny(g1, lb$linked.to, ...))
    o1 && o2
  }
  is.overlap <- sapply(blat.grl, overlapFun, lb=lb, maxgap=500)
  blat.grl <- blat.grl[ is.overlap ]
  blat.grl <- blat.grl [ elementNROWS(blat.grl) == 2 ]
  ##
  ## The split should involve nearly non-overlapping subsequences of
  ## the read, and the total match score should be high (near 100)
  ##
  splitread_ranges <- lapply(blat.grl, sequenceRanges2)
  splitread_int <- splitreadIntersection(splitread_ranges)
  splitread_match <- as.integer(sapply(splitread_ranges, function(x){
    min(sum(x$match), x$Qsize)
  }))
  p <- splitread_match/blat_gr$Qsize[1] * 100
  is_splitread <- splitread_int < 10 & p > 90
  records <- unlist(blat.grl)
  names(records) <- NULL
  ##records_rear <- split(records, records$rear.id)
  ##GRangesList(records_rear)
  records
}



#' Identify rearranged reads -- initiallly unmapped reads that can be
#' aligned by blat to span a novel sequence junction.
#'
#' 
#' For unmapped-mapped read pairs where the mapped read is aligned
#' near a candidate junction, assess whether the unmapped read spans a
#' novel sequence junctions.
#'
#' @return a list of rearranged reads
#' @export
#' @param dirs a \code{DataPaths} object
#' @param rlist a list of \code{RearrangementList} objects
#' @param maxgap a length-one integer vector specifying how much gap
#'   is allowed between a read alignment location and the genomic
#'   intervals for the rearrangement clusters
#' 
rearrangedReadList <- function(dirs, rlist, maxgap=500){
  ids <- names(rlist)
  outfiles <- file.path(dirs[["5parsed_unmapped"]], paste0(ids, ".rds"))
  infiles <- file.path(dirs[["2blat_unmapped"]], paste0(ids, ".rds"))
  results <- setNames(vector("list", length(ids)), ids)
  for(j in seq_along(ids)){
    if(file.exists(outfiles[j])){
      results[[j]] <- readRDS(outfiles[j])
    }
    blat <- readBlat(infiles[j])
    results[[j]] <- rearrangedReads(rlist[[j]], blat, maxgap)
    saveRDS(results[[j]], file=outfiles[j])
  }
  results
}


removeAmbiguous <- function(x, blat){
  blat <- blat[blat$passQC, ]
  x <- x[names(x) %in% blat$rearrangement]
  x
}

#' Removes rearrangements for which the BLAT-aligned reads do not pass QC
#'
#' BLAT records for which the \code{passQC} variable is \code{FALSE}
#' are removed.  The \code{RearrangementList} is then subset to
#' include only those rearrangement ids that remain in the BLAT
#' \code{data.frame}.
#'
#'
#' @seealso See \code{blatScores} for the approach used to QC
#'   BLAT-aligned reads.
#'
#' @export
#' @param rlist a \code{RearrangementList}
#' @param blat_list a list of \code{data.frame} objects, as returned
#'   by \code{scoreBlatExperiment}
#'
#' @examples
#' \dontrun{
#'  library(svovarian)
#'  dirs <- projectOvarian(rootname="OvarianData2")
#'  if(FALSE){
#'  ## A BLAT-filtered RearrangementList
#'  ##saved_result <- readRDS(file.path(dirs[["unit_test"]], "blat_filtered.rds"))
#'  ##tag_list <- readRDS(file.path(dirs[["unit_test"]], "tag_seqs.rds"))
#'  blat <- scoreBlatExperiment(tag_list, dirs)
#'  blat_list <- blat["CGOV2T"]
#'
#'  rlist <- readRDS(file.path(dirs[["rear:filter"]], "CGOV2T.rds"))
#'  rlist <- list(CGOV2T=rlist)
#'  ## Another BLAT-filtered RearrangementList created by 'removeAmbigousAln'
#'  filtered_rlist <- removeAmbiguousAln(rlist, blat_list)
#'  print(filtered_rlist)
#'  }
#' }
removeAmbiguousAln <- function(rlist, blat_list){
  mapply(removeAmbiguous, x=rlist, blat=blat_list)
}