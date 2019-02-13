###############################################################################
# functions in trellis/R/rearrangement-utils.R
###############################################################################
isTranslocation <- function(aln){
  chr1 <- seqnames(first(aln))
  chr2 <- seqnames(last(aln))
  as.logical(chr1!=chr2)
}

R1strand <- function(gpairs) strand(first(gpairs))
R2strand <- function(gpairs) strand(last(gpairs))
R1lessR2 <- function(gpairs){
  is.trans <- isTranslocation(gpairs)
  is.r1less <- as.logical(start(first(gpairs)) < start(last(gpairs)))
  is.r1less[is.trans] <- NA
  is.r1less
}
R1pos <- function(gpairs) as.logical(R1strand(gpairs)=="+")
R2neg <- function(gpairs) as.logical(R2strand(gpairs)=="-")

typeGAPairs <- function(gpairs){
  extdir <- system.file("extdata", package="trellis")
  rear_type <- read.csv(file.path(extdir, "rearrangement_types.csv"), nrows=20,
                        stringsAsFactors=FALSE)
  ##rear_type <- rear_type2
  rtypes <- matrix(NA, length(gpairs), nrow(rear_type))
  is.r1pos <- R1pos(gpairs)
  is.r1neg <- !is.r1pos
  is.r2neg <- R2neg(gpairs)
  is.r2pos <- !is.r2neg
  is.r1less <- R1lessR2(gpairs)
  is.r2less <- !is.r1less
  is.trans <- isTranslocation(gpairs)
  not.trans <- !is.trans
  is.aberrant.sep <- aberrantSep(gpairs)
  not.aberrant.sep <- !is.aberrant.sep
  ##colnames(rtypes) <- rear_type$type
  colnames(rtypes) <- rear_type$name
  rtypes[, 1] <- is.r1pos & is.r2neg & is.r1less & not.trans & not.aberrant.sep
  rtypes[, 2] <- is.r1neg & is.r2pos & is.r2less & not.trans & not.aberrant.sep
  rtypes[, 3] <- is.r1pos & is.r2neg & is.r1less & not.trans & is.aberrant.sep
  rtypes[, 4] <- is.r1neg & is.r2pos & is.r2less & not.trans & is.aberrant.sep
  rtypes[, 5] <- is.r1pos & is.r2neg & is.r1less & not.trans & not.aberrant.sep
  rtypes[, 6] <- is.r1pos & is.r2neg & is.r2less & not.trans
  rtypes[, 7] <- is.r1neg & is.r2pos & is.r2less & not.trans & not.aberrant.sep
  rtypes[, 8] <- is.r1neg & is.r2pos & is.r1less & not.trans
  rtypes[, 9] <-  is.r1pos & is.r2neg & is.trans
  rtypes[, 10] <- is.r1neg & is.r2pos & is.trans
  rtypes[, 11] <- is.r1pos & is.r2neg & is.trans
  rtypes[, 12] <- is.r1neg & is.r2pos & is.trans
  ## inversions
  rtypes[, 13] <- is.r1pos & is.r2pos & is.r1less & not.trans
  rtypes[, 14] <- is.r1pos & is.r2pos & is.r2less & not.trans
  rtypes[, 15] <- is.r1neg & is.r2neg & is.r1less & not.trans
  rtypes[, 16] <- is.r1neg & is.r2neg & is.r2less & not.trans
  rtypes[, 17] <- is.r1pos & is.r2pos &  is.trans
  rtypes[, 18] <- is.r1pos & is.r2pos &  is.trans
  rtypes[, 19] <- is.r1neg & is.r2neg &  is.trans
  rtypes[, 20] <- is.r1neg & is.r2neg &  is.trans
  col.indices <- split(seq_len(ncol(rtypes)), rear_type$compatible)
  y <- vector("list", length(col.indices))
  ##y <- foreach(j = col.indices, .combine="cbind") %do% {
  for(i in seq_along(col.indices)){
    j <- col.indices[[i]]
    y[[i]] <- rowSums(rtypes[, j, drop=FALSE])
  }
  y <- do.call(cbind, y)
  nms <- sapply(split(rear_type$name, rear_type$compatible), function(x) paste(x, collapse=","))
  colnames(y) <- nms
  y
}

rearrangementType <- function(object){
  imp <- improper(object)
  rtypes <- typeGAPairs(imp)
  mns <- pmin(colMeans(rtypes), 1)
  ##nms <- .rearrangement_types()
  nms <- colnames(rtypes)
  data.frame(type=nms[which.max(mns)],
             percent=max(mns), stringsAsFactors=FALSE)
}
###############################################################################
#
###############################################################################
getRearTable <- function() {
  
  type <- c(#paste0(rep("normal", 2), 1:2), 
            paste0(rep("del", 2), 1:2), # deletion
            paste0(rep("amp", 4), 1:4), # amplicon
            paste0(rep("trans", 4), 1:4), # translocation
            paste0(rep("inv", 4), 1:4), # inversion
            paste0(rep("inv", 4), 5:8)) # inverted translocation
  r1 <- c(#"+", "-", # normal
          "+", "-", # deletion
          rep(c("+", "-"), each=2), # amplification
          rep(c("+", "-"), times=2), # translocation
          rep(c("+", "-"), each=2), # inversion
          rep(c("+", "-"), each=2)) # inverted translocation
  r2 <- c(#"-", "+", # normal
          "-", "+", # deletion
          rep(c("-", "+"), each=2), # amplification
          rep(c("-", "+"), times=2), # translocation
          rep(c("+", "-"), each=2), # inversion
          rep(c("+", "-"), each=2)) # inverted translocation
  r1lessr2 <- c(#T, F, # normal
                T, F, # deletion
                T, F, F, T, # amplification
                NA, NA, NA, NA, # translocation
                T, F, T, F, # inversion
                NA, NA, NA, NA) # inverted translocation
  interchrom <- c(#rep(F, 2), # normal
                  rep(F, 6), rep(T, 4), rep(F, 4), rep(T, 4))
  absep <- c(#F, F,
             T, T,
             F, NA, F, NA, 
             rep(NA, 12))
  dnafrag <- c(#"++", "--", 
               "++", "--", 
               "++", "++", "--", "--", 
               "++", "--", "++", "--", 
               "+-", "+-", "-+", "-+", 
               "+-", "+-", "-+", "-+")
  compgroup <- c(#rep(1, 2), 
                 rep(2, 2),
                 rep(c(3, 4), 2), 
                 rep(5, 4),
                 rep(6, 2),
                 rep(7, 2),
                 rep(8,4))
  
  reartypes <- data.frame(type, r1, r2, r1lessr2, interchrom, absep, dnafrag, compgroup)
  return(reartypes)
}

# function that returns junction boundary given rearrangement type
# input: reartype = rearrangement type (e.g. "amp2"), 
#        read = "r1" or "r2",
#        gr = reduced granges object
#        x = margin of interval
# output: list containing chromosome (list$chr) and position (list$pos)
getJunctionBoundary <- function(reartype, read, gr, x) {
  chr <- seqnames(gr)
  
  if(reartype == "del1") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "del2") {pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)}
  else if(reartype == "amp1") {stop('sorry... dunno what to do for amp1')}
  else if(reartype == "amp2") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "amp3") {stop('sorry... dunno what to do for amp3')}
  else if(reartype == "amp4") {pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)}
  else if(reartype == "trans1") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "trans2") {pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)}
  else if(reartype == "trans3") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "trans4") pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)
  else if(reartype == "inv1") {pos <- ifelse(read == "r1", end(gr) + x, end(gr) + x)}
  else if(reartype == "inv2") {pos <- ifelse(read == "r1", end(gr) + x, end(gr) + x)}
  else if(reartype == "inv3") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "inv4") {pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)}
  else if(reartype == "inv5") {pos <- ifelse(read == "r1", end(gr) + x, end(gr) + x)}
  else if(reartype == "inv6") {pos <- ifelse(read == "r1", end(gr) + x, end(gr) + x)}
  else if(reartype == "inv7") {pos <- ifelse(read == "r1", end(gr) + x, start(gr) - x)}
  else if(reartype == "inv8") {pos <- ifelse(read == "r1", start(gr) - x, end(gr) + x)}
  
  return(list("chr" = chr, "pos" = pos))
}


# given galignment pairs of one rearrangement type, return reduced granges of r1s and r2s 
getReducedR1R2 <- function(gap) {
  r1 <- first(gap)
  r2 <- last(gap)
  r1_red <- range(granges(r1))
  r2_red <- range(granges(r2))
  return(list("r1_red" = r1_red, "r2_red" = r2_red))
}

# function to get rearrangement types for given comparable group
# input: comparable group number (numeric)
# output: array of rearrangement types (character)
getTypesFromCompGroup <- function(group) {
  t <- getRearTable()
  types <- t[t$compgroup == group, ]$type
  return(types)
}

# function to get strand of given read 
# GA = GAlignments object
# i = index queried (row/read in GA obj)
getReadStrand <- function(GA, i) as.character(GA[i]@strand@values)

# function that returns opposite strand 
# s = strand, either "+" or "-"
complementStrand <- function(s) {
  if( any(s != "+" & s != "-")) stop('strand must be "+" or "-"')
  return(ifelse(s == "+", "-", "+"))
}

# function that determines if reads in GAlignment object map to different chromosomes
# returns TRUE if interchromosomal, FALSE if not
is.interchrom <- function(GA) length(as.character(GA@seqnames@values)) != 1

# 
# NA matches with both T,F
booleanMatch <- function(b1, b2) {
  if(NA %in% c(b1, b2)) {
    return(TRUE)
  } else {
    return(b1 == b2)
  }
}

# function that determines rearrangement type and comparable group
getType <- function(r1str, r2str, r1lessr2, interchrom, absep) {
  t <- getRearTable()
  rtype <- NA
  compgroup <- NA
  
  # interchrom 
  if(interchrom == TRUE) {
    #print('interchrom = T')
    f <- data.frame(r1str, r2str)
    t1 <- t[t$interchrom==T,]
    # translocation
    if(r1str != r2str) {
      typedf <- t1[grep("trans", t1$type), ]
      rowindex <- prodlim::row.match(f, typedf[, 2:3])
      rtype <- typedf$type[rowindex]
      compgroup <- typedf$compgroup[rowindex]
      
    } else {
      # inverted translocation
      typedf <- t1[grep("inv", t1$type), ]
      rowindex <- prodlim::row.match(f, typedf[, 2:3])
      rtype <- typedf$type[rowindex]
      compgroup <- typedf$compgroup[rowindex]
    }
    
  } else { # same chromosome
    #print("interchrom = F")
    t2 <- t[t$interchrom==F,]
    
    if(r1str == r2str) {
      # inversion
      #print("inversion, r1str == r2str")
      f <- data.frame(r1str, r2str, r1lessr2)
      rowindex <- prodlim::row.match(f, t2[, 2:4])
      rtype <- t2$type[rowindex]
      compgroup <- t2$compgroup[rowindex]
      
    } else { # r1 != r2
      if(absep == T) {
        #print("deletion1,2 or amp2,4")
        # deletion(1,2) or amp(2,4)
        t2a <- t2[c(1,2,4,6), ]
        f <- data.frame(r1str, r2str, r1lessr2, interchrom)
        rowindex <- prodlim::row.match(f, t2a[, 2:5])
        rtype <- t2a$type[rowindex]
        compgroup <- t2a$compgroup[rowindex]
        
      } else {
        #print("amp1,2,3,4")
        # amp(1,2,3,4)
        t2b <- t2[grep("amp", t2$type), ]
        f <- data.frame(r1str, r2str, r1lessr2, interchrom)
        rowindex <- prodlim::row.match(f, t2b[, 2:5])
        rtype <- t2b$type[rowindex]
        compgroup <- t2b$compgroup[rowindex]
      }
    }
  }
  return(list("rtype" = rtype, "compgroup" = compgroup))
}

# function to summarize improper read pairs (GAlignmentPairs object)
# input: GAlignmentPairs object
# output: GAlignmentPairs object with mcols:
#   r1_strand   r2_strand  r1lessr2 dna_fragment interchromosomal aberrant_separation rearrangement_type comparable_group
summarizeImproperRP <- function(irp) {
  
  irp_withmeta <- irp
  
  r1str_all <- as.character(R1strand(irp))
  r2str_all <- as.character(R2strand(irp))
  dnafrag_all <- paste0(r1str_all, complementStrand(r2str_all))
  interchrom_all <- unname(sapply(irp, is.interchrom))
  ab_sep_all <- aberrantSep(irp, distance=10e3)
  r1lessr2_all <- R1lessR2(irp)
  
  # map to rearrangement type/comparable group
  t <- getRearTable()
  rtype <- c()
  compgroup <- c()
  for(i in 1:length(irp)) {
    info <- getType(r1str_all[i], r2str_all[i], r1lessr2_all[i], interchrom_all[i], ab_sep_all[i])
    rtype <- c(rtype, info$rtype)
    compgroup <- c(compgroup, info$compgroup)
  }
  
  mcols(irp_withmeta) <- data.frame(r1_strand = r1str_all, 
                                    r2_strand = r2str_all,
                                    r1lessr2 = r1lessr2_all,
                                    dna_fragment = dnafrag_all,
                                    interchromosomal = interchrom_all,
                                    aberrant_separation = ab_sep_all,
                                    rearrangement_type = rtype,
                                    comparable_group = compgroup)
  return(irp_withmeta)
}

getCompGroupFromSIRP <- function(sirp) {
  cg_array <- sirp@elementMetadata$comparable_group
  cg_tab <- table(cg_array)
  return(names(cg_tab)[which.max(cg_tab)])
}

# take 2, reduce first then label
summarizeIRP <- function(irp) {
  
  irp_withmeta <- irp
  
  r1str_all <- as.character(R1strand(irp))
  r2str_all <- as.character(R2strand(irp))
  dnafrag_all <- paste0(r1str_all, complementStrand(r2str_all))
  interchrom_all <- unname(sapply(irp, is.interchrom))
  ab_sep_all <- aberrantSep(irp, distance=10e3)
  r1lessr2_all <- R1lessR2(irp)
  
  # map to rearrangement type/comparable group
  t <- getRearTable()
  rtype <- c()
  compgroup <- c()
  for(i in 1:length(irp)) {
    info <- getType(r1str_all[i], r2str_all[i], r1lessr2_all[i], interchrom_all[i], ab_sep_all[i])
    rtype <- c(rtype, info$rtype)
    compgroup <- c(compgroup, info$compgroup)
  }
  
  mcols(irp_withmeta) <- data.frame(r1_strand = r1str_all, 
                                    r2_strand = r2str_all,
                                    r1lessr2 = r1lessr2_all,
                                    dna_fragment = dnafrag_all,
                                    interchromosomal = interchrom_all,
                                    aberrant_separation = ab_sep_all,
                                    rearrangement_type = rtype,
                                    comparable_group = compgroup)
  return(irp_withmeta)
}


