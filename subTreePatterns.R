# calculate presence absence data for subtrees
load("~/data/IES_data/rdb/charMats")
source("~/projects/IES/src/sharedFunctions.R")
sp0 <- read.table("~/data/IES_data/rdb/Pca.PbiPtePse.dat", as.is = TRUE, header = TRUE)
sp2 <- read.table("~/data/IES_data/rdb/Pse.PbiPte.dat", as.is = TRUE, header = TRUE)
patternsV0 <- presenceAbsence(sp0)
patternsV2 <- presenceAbsence(sp2)
save(patternsV2, patternsV0, file = "~/data/IES_data/rdb/subTreePatterns")

