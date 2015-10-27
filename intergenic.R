# find intergenic homologs 

source("~/projects/IES/src/sharedFunctions.R")

# load all pairs of genes
load("~/data/IES_data/intergenic/genePairs")
homgroups <- read.table("~/data/IES_data/rdb/homologGroups.dat", stringsAsFactors = FALSE, header = TRUE)

# keep only genes in homology groups
pbi <- pbiGP[(which(pbiGP[, 1] %in% homgroups$geneId)), ]
pte <- pteGP[(which(pteGP[, 1] %in% homgroups$geneId)), ]
pse <- pseGP[(which(pseGP[, 1] %in% homgroups$geneId)), ]
pca <- pcaGP[(which(pcaGP[, 1] %in% homgroups$geneId)), ]

allSpecM <- rbind(pbi, pte, pse, pca)
homgroupOfGenesI <- match(allSpecM[, 1], homgroups[, "geneId"])
homgroupOfNeighborsI <- match(allSpecM[, 2], homgroups[, "geneId"])
extM <- cbind(allSpecM[, 1:5], homgroups[homgroupOfGenesI, "id"], homgroups[homgroupOfNeighborsI, "id"], allSpecM[, 6:8])
hgV <- unique(extM[, 6]) # vector of homologous groups
synblockId <- 1 # id of syntenic blocks
intergenicDF <- data.frame("synblockId"   = numeric(),
                         "upanchor"     = character(),
                         "downanchor"   = character(),
                         "length"       = numeric(),
                         "uphomgroup"   = numeric(),
                         "downhomgroup" = numeric(),
                         "begin"        = numeric(),
                         "end"          = numeric(),
                         "scaffold"     = character())
counter <- 1
noChoice <- 1 # counter for the number of times a homolog cannot be chosen (from two alternatives)
for(hg in hgV){
  #hg <- hgV[1]
  cat(counter,"/", length(hgV), "\r")
  # ignore genes with neighbors not in homologous groups (NA)
  naIndex <- is.na(extM[, 7])
  hgI <- intersect(which(extM[, 6] == hg), which(!naIndex))
  if(length(hgI) == 1){ # not a block if there is only one pair of anchors
    next
  }
  # make a key of rel. orientation and group
  oriGroup <- apply(extM[hgI, c(5,7)], 1, paste, collapse = " ") 
  # find if there are more than one neighbors that belnong to the same homology group and have compatible orientations
  allduplicated <- duplicated(oriGroup) | duplicated(oriGroup, fromLast = TRUE)
  # homologous groups with more than one links
  hgWithLinks <- unique(extM[hgI[allduplicated], 7])
  if(any(duplicated(extM[hgI[allduplicated], 2]))){
    noChoice <- noChoice + 1
    next
  }
  block <- extM[hgI[allduplicated], ]
  for(synblock in hgWithLinks){
    #synblock <- hgWithLinks[1]
    DF <- data.frame(synblock, block[block[, 7] == synblock, c(1, 2, 4, 6, 7, 8, 9, 10)], stringsAsFactors = FALSE)
    names(DF) <- c("synblockId", "upanchor", "downanchor", "length", "uphomgroup", "downhomgroup", "begin", "end", "scaffold")
    intergenicDF <- rbind(intergenicDF, DF)
  }
  synblockId <- synblockId + 1
  counter <- counter + 1
}
cat(paste("no choice made over homology in", noChoice, "cases\n"))
save(intergenicDF, noChoice, file = "~/data/IES_data/rdb/intergenic")
