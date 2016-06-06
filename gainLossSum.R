# adjust gain and loss rate by size of alignment and summarize resutls

gl <- read.table("~/data/IES/analysis/tables/gainLoss.dat", header = TRUE)

# calculate total length of conserved blocks of alignments
gb <- read.table("~/data/IES/analysis/tables/gblocks.dat")
geneFamilies <- unique(gb$V1)
blockLenghts <- gb$V3 - gb$V2 + 1
bL <- data.frame(geneFamily = geneFamilies, nucleotides = numeric(length(geneFamilies)))
for(geneFamily in geneFamilies){
  gfI <- which(gb$V1 == geneFamily)
  bL$nucleotides[which(bL$geneFamily == geneFamily)] <- sum(blockLenghts[gfI])
}

key <- paste(gl$from, gl$to, sep = '-')
Ploss <- gl$Panc * gl$loss
nt <- bL$nucleotides[match(gl$cluster, bL$geneFamily)]

l <- length(unique(key))
DF <- data.frame(branch = character(l), gainRate = numeric(l), lossRate = numeric(l), stringsAsFactors = FALSE)
counter <- 1
for(k in unique(key)){
  index <- which(key == k)
  pgain <- sum(gl$gain[index]/nt[index])/length(index)
  #pgain <- sum(gl$gain[index])/length(index)
  ploss <- sum(Ploss[index])/length(index)
  DF[counter, 1] <- k
  DF[counter, 2:3] <- c(pgain, ploss)
  counter <- counter + 1
}
#save(DF, file = "~/data/IES_data/rdb/ratesPerBranch")
DF$gainRate <- DF$gainRate * 10^5
DF$lossRate <- DF$lossRate * 10^2
knitr::kable(DF[order(DF$branch), ], row.names = FALSE, digits = 2)
