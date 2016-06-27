#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# calculate the normalize branch lengths for single gene family trees

if (length(args) != 2) {
  stop("Missing argument(s), input output")
}
geneTreeSpeciesTreeF <- args[1] #~/data/IES/analysis/tables/geneTreeSpeciesTree.tab"
normBrLensF <- args[2] #"~/data/IES/analysis/tables/normBrLens.tab"

library(ape)
source("~/projects/IES/src/sharedFunctions.R")

nodeDictionary <- read.table("~/data/IES/analysis/tables/nodeDictionary.dat", stringsAsFactors = FALSE, header = TRUE)
spEvent <- read.table("~/data/IES/analysis/tables/spEvents.dat", header = TRUE, stringsAsFactors = FALSE)
singleGeneFamilies <- unname(unlist(read.table(geneTreeSpeciesTreeF, header = FALSE, stringsAsFactors = FALSE)))

# branches of interest
boe <- c("b0_1",
         "b0_2",
         "b1_4",
         "b4_5",
         "b5_7",
         "b5_8",
         "b4_6",
         "b6_9",
         "b9_11",
         "b11_13",
         "b11_14",
         "b9_12",
         "b12_16",
         "b16_17",
         "b16_18",
         "b12_15",
         "b6_10",
         "b1_3")
m <- matrix(ncol = length(boe), nrow = length(singleGeneFamilies))
sgfs <- character() #selected gene families
counter <- 1
c1 <- 1
for(gf in singleGeneFamilies){
  # read tree
  tr <- read.tree(paste("~/data/IES/analysis/phyldog/results/", gf, ".ReconciledTree", sep = ""))
  # proceed only if gene family has IES
  c1 <- c1 + 1
  gfI <- which(spEvent$cluster == gf)
  if(length(gfI) == 0){
    # was not included in the ancestral state reconstruction, e.g. does not have IESs
    next
  }
  print(paste(counter, "/", c1, "/", length(singleGeneFamilies)))
  totalBrL <- sum(tr$edge.length) # normalize branch length by total length
  DF <- data.frame(tr$edge, nbrl = tr$edge.length/totalBrL)
  # translate node id from r to phyldog
  DF$X1 <- r2phyldog(DF$X1, gf)
  DF$X2 <- r2phyldog(DF$X2, gf)
  nd2e <- spEvent[gfI, ]
  DF$X1 <- nd2e$spEvent[match(DF$X1, nd2e$nodeP)]
  DF$X2 <- nd2e$spEvent[match(DF$X2, nd2e$nodeP)]
  key <-  paste("b", DF$X1, "_",  DF$X2, sep = "")
  m[counter, ] <- DF[match(boe, key), "nbrl"]
  sgfs <- append(sgfs, gf)
  counter <- counter + 1
}
m <- m[rowSums(is.na(m)) != ncol(m),]

d <- data.frame(m, row.names = sgfs, stringsAsFactors = FALSE)
# keep gene families that are within the 95th percentile of all branghes
names(d) <- boe
write.table(d, file = normBrLensF, sep = "\t")
