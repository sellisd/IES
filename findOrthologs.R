# find pairs of orthologs
# read pairs of true orthologs
orth <- read.table("~/data/IES_data/msas/orthologs/pbipte.dat", as.is = TRUE, header=TRUE)
# make a convenience key from both names pasted
key <- paste(orth$pbi, orth$pte)
# read character matrices
path <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(path, pattern="cluster.[^.]*.dat")
matBool <- matrix(nrow = windowSize * 2, ncol = 0)
logoList <- list()
logoInfo <- character()
for(fileName in files){
  #  fileName <- files[7]
  a <- read.table(paste("~/data/IES_data/msas/alignments/charMat/", fileName,sep=""), header=T, row.names=1, as.is=TRUE)
  # find which species are represented
  speciesV <- substr(row.names(a), 0,4)
  pbiI <- which(speciesV == "PBIA")
  pteI <- which(speciesV == "PTET")
  if(length(pbiI)==0 | length(pteI)==0){
    next # skip alignments without both P biaurelia and P. tetraurelia genes
  }
  # find rows that IES sequence has to be reverse complemented
  toRevComp <- c(which(pbiStrands[row.names(a)] == "-"), which(pteStrands[row.names(a)] == "-"))
  for(columnI in c(1:ncol(a))){ # for each column
    #columnI <- 1
    allInColumn <- character() # all homologous strings
    presentI <- which(a[, columnI]!=0)
    pbiRows <- intersect(presentI, pbiI)
    pteRows <- intersect(presentI, pteI)
    if(length(pbiRows) == 0 | length(pteRows) == 0){
      next # skip columns without both P. biaurelia and P. tetraurelia IES
    }
    queryPairs <- combS(pbiRows, pteRows) # make all combinations of pbiRows and pteRows
    trueOrthologsI <- queryPairs[paste(row.names(a)[queryPairs[, 1]],row.names(a)[queryPairs[, 2]]) %in% key, , drop = FALSE]
    orthIES <- data.frame(pbi = a[trueOrthologsI[,1], columnI], pte = a[trueOrthologsI[,2], columnI]) # pairs of orthologous IES
    # >>> Here we could possibly keep randomly only one of each IES (remove random duplicates)
    discardRow <- pteD$floating[orthIES$pte] | pbiD$floating[orthIES$pbi] # if one is floating discard both
    orthIES <- orthIES[!discardRow, ]
    trueOrthologsI <- trueOrthologsI[which(!discardRow), , drop = FALSE]
    if(nrow(orthIES) == 0){
      next  # skip if no pairs are left
    }
    # >>> Here check 10bp upstream and downstream there is only CD
    # for each row find if they need reverse complement and pull front , back and flanking sequences
    for(i in c(1:nrow(trueOrthologsI))){
      pbiFront <- pbiD[orthIES[i, 1], "front"]
      pbiBack <- pbiD[orthIES[i, 1], "back"]
      pteFront <- pteD[orthIES[i, 2], "front"]
      pteBack <- pteD[orthIES[i, 2], "back"]
      if(trueOrthologsI[i, 1] %in% toRevComp){
        pbiMV <- rev(comp(s2c(paste0(pbiFront,pbiBack)))) # swap front back and revecse complement
      }else{
        pbiMV <- s2c(paste0(pbiFront,pbiBack))
      }
      if(trueOrthologsI[i, 2] %in% toRevComp){
        pteMV <- rev(comp(s2c(paste0(pteFront,pteBack))))
      }else{
        pteMV <- s2c(paste0(pteFront,pteBack))
      }
      matBool <- cbind(matBool, pairwiseIdentity(pbiMV, pteMV))
      # add to a character matrix and make a sequenceLogo for each pair
      allInColumn <- append(allInColumn, toupper(c(c2s(pbiMV), c2s(pteMV))))
    }
    logoList[length(logoList)+1] <- list(allInColumn)
    logoInfo <- append(logoInfo, paste(fileName))
    #prepareLOGO(toupper(allInColumn))
  }
  #  aln <- as.alignment(nb=length(allInColumn), seq = toupper(allInColumn))
  #if(length(allInColumn) == 0){print "haho"; stop()}
}
