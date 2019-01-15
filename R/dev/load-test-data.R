library(Rsamtools)
library(magrittr)
library(Rsamtools)
library(ggplot2)
library(svbams)
library(svfilters.hg19)
library(tidyverse)
library(gridExtra)
library(trellis)
library(ggnet)
library(graph)
library(Rgraphviz)
library(network)
library(sna)
library(kableExtra)
library(GenomicAlignments)

# load data from mini bam
extdata <- system.file("extdata", package="svbams")
bamfile <- file.path(extdata, "cgov44t_revised.bam")
what <- c("flag", "mrnm", "mpos", "mapq")
iparams <- improperAlignmentParams(what=what)
improper_rp <- getImproperAlignmentPairs(bamfile,
                                         param=iparams,
                                         build="hg19")
data(bins1kb, package="svfilters.hg19")
bins <- keepSeqlevels(bins1kb, "chr15",
                      pruning.mode="coarse")
flags <- scanBamFlag(isDuplicate=FALSE,
                     isPaired=TRUE,
                     isUnmappedQuery=FALSE,
                     hasUnmappedMate=FALSE,
                     isProperPair=TRUE)
bviews <- BamViews(bamPaths=bamfile, bamRanges=bins)
bins$cnt <- binnedCounts(bviews)
bins <- bins[ bins$cnt > 0 ]
bins$std_cnt <- binNormalize(bins)
set.seed(123)
bins$log_ratio <- binGCCorrect(bins)
g <- segmentBins(bins, param=SegmentParam())
g
data(bins1kb, package="svfilters.hg19")
ddir <- system.file("extdata", package="svbams",
                    mustWork=TRUE)
lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
bins1kb$log_ratio <- lr
bins <- keepSeqlevels(bins1kb, c("chr5", "chr8", "chr15"),
                      pruning.mode="coarse")
path <- system.file("extdata", package="svbams")
segs <- readRDS(file.path(path, "cgov44t_segments.rds"))
seqlevels(segs, pruning.mode="coarse") <- seqlevels(bins)

# deletions
dp <- DeletionParam(remove_hemizygous=FALSE)
dp
del.gr <- IRanges::reduce(segs[segs$seg.mean < hemizygousThr(dp)],
                          min.gapwidth=2000)
proper_rp <- properReadPairs(bamfile, gr=del.gr, dp)
improper_rp <- keepSeqlevels(improper_rp, seqlevels(segs),
                             pruning.mode="coarse")
read_pairs <- list(proper_del=proper_rp, improper=improper_rp)
pdata <- preprocessData(bam.file=bamfile,
                        genome="hg19",
                        bins=bins,
                        segments=segs,
                        read_pairs=read_pairs)
deletions <- sv_deletions(preprocess=pdata, param=dp)
variant(deletions)
calls(deletions)
improper(deletions[[4]])

# rearrangements
rparam <- RearrangementParams(min_number_tags_per_cluster=5,
                              rp_separation=10e3)
minNumberTagsPerCluster(rparam)
rpSeparation(rparam) %>% "/"(1000) %>% paste0("kb")
rlist <- findCandidates2(pdata, rparam)
rlist