context("Split reads")

test_that("split_reads", {
  data(rlist, package="trellis")
  extdata <- system.file("extdata", package="svbams")
  unmap.file <- file.path(extdata, "blat_unmapped.txt")
  blat_unmap <- readBlat(unmap.file)
  split_reads <- rearrangedReads(linkedBins(rlist), blat_unmap, 500)
  expect_identical(names(split_reads), names(rlist))
  expect_is(split_reads, "GRangesList")
  splitReads(rlist) <- split_reads
  expect_identical(split_reads, splitReads(rlist))
})
