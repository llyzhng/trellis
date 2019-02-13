
calc_VAF <- function(normalReads, tumorReads){
  return (tumorReads / (normalReads/2 + tumorReads))
}

# function to count number of normal reads in the interval [position +/- x] in given bam file
# input: bamfile, genomic position, x=margin size
# output: number of normal reads 
countNormalReads <- function(bamfile, chr, pos, x) {
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = pos - x, 
                                 end = pos + x))
  nparam <- ScanBamParam(which = gr, 
                         flag = scanBamFlag(isProperPair = T, 
                                            isPaired = T))
  cb <- countBam(bamfile, param=nparam)
  return(cb$records)
}

# function to annotate deletions granges with counts of normal reads
# output granges with normal_reads mcol
annotateDeletionsNormal <- function(bamfile, del_gr, x=100) {
  #del_gr <- variant(deletions)
  
  # count normal reads at junctions
  normal_reads <- c()
  for (i in 1:length(del_gr)) {
    # left boundary
    left_reads <- countNormalReads(bamfile = bamfile, 
                                   chr = seqnames(del_gr[i]),
                                   pos = start(del_gr[i]),
                                   x = x)
    # right boundary interval
    right_reads <- countNormalReads(bamfile = bamfile, 
                                   chr = seqnames(del_gr[i]),
                                   pos = end(del_gr[i]),
                                   x = x)

    normal_reads <- c(normal_reads, sum(left_reads, right_reads))
    
    # improper read pairs
    del_irp <- deletions[[i]]
    #summarizeImproperRP(del_irp)
  }
  mcols(del_gr)$normal_reads <- normal_reads
  return(del_gr)
}



###############################################################################
###############################################################################

# function to calculate VAF of SV
# input: Rearrangement object (e.g. rlist[[1]])
calcVAFFromIRP <- function(bamfile, x=100, irp) {
  sum <- summarizeImproperRP(irp)
  
  # which comparable group
  group <- names(table(sum@elementMetadata$comparable_group))[which.max(table(sum@elementMetadata$comparable_group))]
  
  # types of rearrangements for that comparable group
  types <- getTypesFromCompGroup(as.numeric(group))
  
  tumor_count <- c()
  normal_count <- c()

  for(t in types) {
    gap <- sum[sum@elementMetadata@listData$rearrangement_type == t]
    tcount <- length(gap)
    red <- getReducedR1R2(gap)
    junc_r1 <- getJunctionBoundary(reartype = t, read = "r1", gr = red$r1_red, x = x)
    junc_r2 <- getJunctionBoundary(reartype = t, read = "r2", gr = red$r2_red, x = x)
    n1 <- countNormalReads(bamfile, chr=junc_r1$chr, junc_r1$pos, x)
    n2 <- countNormalReads(bamfile, chr=junc_r2$chr, junc_r2$pos, x)
    ncount <- sum(n1, n2)
    
    tumor_count <- c(tumor_count, tcount)
    normal_count <- c(normal_count, ncount)
  }
  
  VAF <- calc_VAF(sum(normal_count), sum(tumor_count))
  return(VAF)
}

calcVAFFromIRP_del <- function(bamfile, x=100, irp) {
  sum <- summarizeImproperRP(irp)
  
  # which comparable group
  group <- names(table(sum@elementMetadata$comparable_group))[which.max(table(sum@elementMetadata$comparable_group))]
  # correct if mislabeled as amp1,3, compgroup3
  if(group == 3) {
    sum@elementMetadata$aberrant_separation <- TRUE
    rtype <- sum@elementMetadata$rearrangement_type
    sum@elementMetadata$rearrangement_type <- ifelse(rtype == "amp1", "del1", ifelse(rtype == "amp3", "del2", NA))
    sum@elementMetadata$comparable_group <- 2
  } else if(group == 2) {
    continue
  } else {
    stop(paste0("something is wrong...compgroup=",group))
  }
  
  tumor_count <- c()
  normal_count <- c()
  
  for(t in c("del1", "del2")) {
    gap <- sum[sum@elementMetadata@listData$rearrangement_type == t]
    tcount <- length(gap)
    red <- getReducedR1R2(gap)
    junc_r1 <- getJunctionBoundary(reartype = t, read = "r1", gr = red$r1_red, x = x)
    junc_r2 <- getJunctionBoundary(reartype = t, read = "r2", gr = red$r2_red, x = x)
    n1 <- countNormalReads(bamfile, chr=junc_r1$chr, junc_r1$pos, x)
    n2 <- countNormalReads(bamfile, chr=junc_r2$chr, junc_r2$pos, x)
    ncount <- sum(n1, n2)
    
    tumor_count <- c(tumor_count, tcount)
    normal_count <- c(normal_count, ncount)
  }
  
  VAF <- calc_VAF(sum(normal_count), sum(tumor_count))
  return(VAF)
}

calc_VAF_amp24<- function(bamfile, x=100, Rearrangement) {
  irp <- improper(Rearrangement)
  sum <- summarizeImproperRP(irp)
  # comp group 4: amp2,4
  
  # amp2
  amp2_gap <- sum[sum@elementMetadata@listData$rearrangement_type == "amp2"]
  amp2_red <- getReducedR1R2(amp2_gap)
  amp2_junc_r1 <- getJunctionBoundary(reartype = "amp2", read = "r1", gr = amp2_red$r1_red, x = x)
  amp2_junc_r2 <- getJunctionBoundary(reartype = "amp2", read = "r2", gr = amp2_red$r2_red, x = x)
  
  amp2_n1 <- countNormalReads(bamfile, chr=amp2_junc_r1$chr, amp2_junc_r1$pos, x)
  amp2_n2 <- countNormalReads(bamfile, chr=amp2_junc_r2$chr, amp2_junc_r2$pos, x)
  amp2_t <- length(amp2_gap)
  amp2_n <- sum(amp2_n1, amp2_n2)
  
  # amp4
  amp4_gap <- sum[sum@elementMetadata@listData$rearrangement_type == "amp4"]
  amp4_red <- getReducedR1R2(amp4_gap)
  amp4_junc_r1 <- getJunctionBoundary(reartype = "amp4", read = "r1", gr = amp4_red$r1_red, x = x)
  amp4_junc_r2 <- getJunctionBoundary(reartype = "amp4", read = "r2", gr = amp4_red$r2_red, x = x)
  
  amp4_n1 <- countNormalReads(bamfile, chr=amp4_junc_r1$chr, amp4_junc_r1$pos, x)
  amp4_n2 <- countNormalReads(bamfile, chr=amp4_junc_r2$chr, amp4_junc_r2$pos, x)
  amp4_t <- length(amp4_gap)
  amp4_n <- sum(amp4_n1, amp4_n2)
  
  VAF <- calc_VAF(sum(amp2_n, amp4_n), sum(amp2_t, amp4_t))
  return(VAF)
}