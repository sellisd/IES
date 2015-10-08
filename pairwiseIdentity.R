# calculate pairwise identities for each position in orthologous IES between P. biaurelia and P. tetraurelia

# set default value
windowSize <- 20

# load necessary data
load("~/data/IES_data/rdb/IESorthologs") 
load("~/data/IES_data/rdb/conservationPatterns")
load("~/data/IES_data/rdb/lengthBins")
source("~/projects/IES/src/sharedFunctions.R")

# prepare boolean matrices to fill in with the results of pairwise identities
matBool <- matrix(nrow = windowSize * 2, ncol = 0)
matBoolR <- matrix(nrow = windowSize * 2, ncol = 0) # backround randomized
# extend pairs of orthologs table to include length bin and conservation pattern
key <- paste(cons$group, cons$col)
extendedTable <- cbind(orthPairs, lengthBin = pteLengthBins[orthPairs$pte], pattern = cons$pattern[match(orthPairs$cluster, key)])
# permute P. tetraurelia IES keeping fixed the length bin assignment AND conservation pattern when available
permutedI <- permuteBy(extendedTable, c("pte", "pteseq"), c("lengthBin", "pattern"),ret = "index")
# There are cases where we do not have a 1-1 orthology assignment, treat clusters all together 
clusters <- unique(orthPairs$cluster) #  cluster column
# the IES sequence is 2bp larger than the actual length as it contains both flanking TAs
pbiL <- nchar(orthPairs$pbiseq)
pteL <- nchar(orthPairs$pteseq)
pteRL <- nchar(orthPairs$pteseq[permutedI])
# split out front and back sequence
frontPbi <- substr(orthPairs$pbiseq, 0, windowSize)
frontPte <- substr(orthPairs$pteseq, 0, windowSize)
frontPteR <- substr(orthPairs$pteseq[permutedI], 0, windowSize)
backPbi  <- substr(orthPairs$pbiseq, pbiL - windowSize + 1, pbiL)
backPte  <- substr(orthPairs$pteseq, pteL - windowSize + 1, pteL)
backPteR <- substr(orthPairs$pteseq[permutedI], pteRL - windowSize + 1, pteRL)

for(i in clusters){
  #i <- clusters[1]
  index <- which(orthPairs$cluster == i)
  if(length(index) == 0){
    next
  }
  borderSeqPbi <- paste0(frontPbi[index], backPbi[index])
  borderSeqPte <- paste0(frontPte[index], backPte[index])
  borderSeqPteR <- paste0(frontPteR[index], backPteR[index])
  for(j in c(1:length(borderSeqPbi))){
    # calculate pairwise sequence identity
    matBool <- cbind(matBool, pairwiseIdentity(s2c(borderSeqPbi[j]), s2c(borderSeqPte[j])))
    matBoolR <- cbind(matBoolR, pairwiseIdentity(s2c(borderSeqPbi[j]), s2c(borderSeqPteR[j])))
  }
}

save(matBool, matBoolR, extendedTable, file = "~/data/IES_data/rdb/pairwiseIdentity")
