# find pairs of orthologs
library( seqinr)
load("~/data/IES_data/rdb/iesInfo")
load("~/data/IES_data/rdb/geneStrand")
source("~/projects/IES/src/sharedFunctions.R")
# read pairs of true orthologs
orth <- read.table("~/data/IES_data/msas/orthologs/pbipte.dat", as.is = TRUE, header=TRUE)
orthPairs <- data.frame()
# make a convenience key from both names pasted
key <- paste(orth$pbi, orth$pte)
# read character matrices
path <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(path, pattern="cluster.[^.]*.dat")
for(fileName in files){
  #  fileName <- files[7]
  # fileName <- "cluster.5591.dat"
  #fileName <- "cluster.10008.dat"
  a <- read.table(paste(path, fileName,sep=""), header=T, row.names=1, as.is=TRUE)
  # find which species are represented
  speciesV <- substr(row.names(a), 0,4)
  pbiI <- which(speciesV == "PBIA")
  pteI <- which(speciesV == "PTET")
  if(length(pbiI)==0 | length(pteI)==0){
    next # skip alignments without both P biaurelia and P. tetraurelia genes
  }
  # find rows that IES sequence has to be reverse complemented
  toRevComp <- c(which(pbiStrand[row.names(a)] == "-"), which(pteStrand[row.names(a)] == "-"))
  for(columnI in c(1:ncol(a))){ # for each column
    #columnI <- 2
    allInColumn <- character() # all homologous strings
    presentI <- which(a[, columnI]!=0)
    pbiRows <- intersect(presentI, pbiI)
    pteRows <- intersect(presentI, pteI)
    if(length(pbiRows) == 0 | length(pteRows) == 0){
      next # skip columns without both P. biaurelia and P. tetraurelia IES
    }
    queryPairs <- combS(pbiRows, pteRows) # make all combinations of pbiRows and pteRows
    trueOrthologsI <- queryPairs[paste(row.names(a)[queryPairs[, 1]],row.names(a)[queryPairs[, 2]]) %in% key, , drop = FALSE]
    orthIES <- data.frame(pbi = a[trueOrthologsI[,1], columnI], pte = a[trueOrthologsI[,2], columnI], stringsAsFactors = FALSE) # pairs of orthologous IES
    # >>> Here we could possibly keep randomly only one of each IES (remove random duplicates)
    discardRow <- pteD[orthIES$pte, "floating"] | pbiD[orthIES$pbi, "floating"] # if one is floating discard both
    orthIES <- orthIES[!discardRow, , drop = FALSE]
    trueOrthologsI <- trueOrthologsI[which(!discardRow), , drop = FALSE]
    if(nrow(orthIES) == 0){
      next  # skip if no pairs are left
    }
    cluster <- as.integer(sub(".dat","",sub("cluster.","",fileName)))
    for(i in c(1:nrow(trueOrthologsI))){
      #i <- 1
      pbiIESSeq <- s2c(pbiD[orthIES[i,1], "seq"])
      pteIESSeq <- s2c(pteD[orthIES[i,2], "seq"])
      if(trueOrthologsI[i, 1] %in% toRevComp){
        pbiIESSeq <- rev(comp(pbiIESSeq))
      }
      if(trueOrthologsI[i, 2] %in% toRevComp){
        pteIESSeq <- rev(comp(pteIESSeq))
      }
      if(is.na(orthIES[i,1])){
        print(cluster)
        next()
        # skip the rare cases of two IES on the same column (not sure if it is one or two)
        # currently the only known example is of 5591
      }
      orthPairs <- rbind(orthPairs, data.frame(cluster = paste(cluster, columnI), pbi = orthIES[i,1], pte = orthIES[i,2], pbiseq = toupper(c2s(pbiIESSeq)), pteseq = toupper(c2s(pteIESSeq)), stringsAsFactors = FALSE))
    }
  }
}
save(orthPairs, file = "~/data/IES_data/rdb/IESorthologs")
