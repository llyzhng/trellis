suppressMessages(source("~/GitHub/trellis/R/dev/load-test-data.R"))
source("~/GitHub/trellis/R/dev/rear-functions.R")
source("~/GitHub/trellis/R/dev/calc-SV-VAF.R")

# interval around breakpoint (junction +/- x)
x <- 100

###############################################################################
# rearrangements
rlist
rlist[[1]]
irp1 <- improper(rlist[[1]])
irp1

sum <- summarizeImproperRP(irp1)
sum

calc_VAF_amp24(bamfile, x, rlist[[1]])
calc_VAF_amp24(bamfile, x, rlist[[2]])

calcVAFFromIRP(bamfile, x, irp1)
calcVAFFromIRP(bamfile, x, improper(rlist[[2]]))

###############################################################################
# deletions data 

deletions
variant(deletions)
calls(deletions)
del_irp <- improper(deletions[[4]])

calcVAFFromIRP_del(bamfile, x, del_irp)

sum <- summarizeImproperRP(del_irp)
sum
t <- getRearTable()

#calcVAFFromIRP(bamfile, x, del_irp)
#annotateDeletionsNormal(bamfile, variant(deletions), x=100)


