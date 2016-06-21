# adjust gain and loss rate by size of alignment and summarize resutls
source("~/projects/IES/src/sharedFunctions.R")

gainloss <- read.table("~/data/IES/analysis/tables/gainLoss.dat", header = TRUE)

# calculate total length of conserved blocks of alignments
gb <- read.table("~/data/IES/analysis/tables/gblocks.dat")
geneFamilies <- unique(gb$V1)
blockLenghts <- gb$V3 - gb$V2 + 1
bL <- data.frame(geneFamily = geneFamilies, nucleotides = numeric(length(geneFamilies)))
for(geneFamily in geneFamilies){
  gfI <- which(gb$V1 == geneFamily)
  bL$nucleotides[which(bL$geneFamily == geneFamily)] <- sum(blockLenghts[gfI])
}

key <- paste(gainloss$from, gainloss$to, sep = '-')
Ploss <- gainloss$Panc * gainloss$loss
nt <- bL$nucleotides[match(gainloss$cluster, bL$geneFamily)]

l <- length(unique(key))
DF <- data.frame(fromS = character(l), toS = character(l), gainRate = numeric(l), lossRate = numeric(l), stringsAsFactors = FALSE)
counter <- 1
for(k in unique(key)){
  index <- which(key == k)
  pgain <- sum(gainloss$gain[index]/nt[index])/length(index)
  #pgain <- sum(gainloss$gain[index])/length(index)
  ploss <- sum(Ploss[index])/length(index)
  DF[counter, 1] <- gainloss$from[index][1]
  DF[counter, 2] <- gainloss$to[index][1]
  DF[counter, 3:4] <- c(pgain, ploss)
  counter <- counter + 1
}

# adjust by branch length manually for now!
brlt <- read.tree("~/data/IES/analysis/brlen/part.nexus.treefile") # all models where CODON models, branch lengths are number of substitutions per codon site!
nodelabels()
tiplabels()
brlensDF <- data.frame(brlt$edge, brlt$edge.length / 3, stringsAsFactors = FALSE) # rescaled to number of substitutions per nt site

rDF <- DF # round and scale for plotting
rDF$gainRate <- paste("+", round(DF$gainRate * 10^5, 2))
rDF$lossRate <- paste("-", round(DF$lossRate * 10^2, 2))

write.table(rDF, file = "~/data/IES/analysis/tables/ratesPerBranch.dat", row.names = FALSE, quote = FALSE, sep = "\t")

knitr::kable(DF[order(DF$branch), ], row.names = FALSE, digits = 2)
