# read probability of presence in each gene tree nodepath (in revBayes nodeId)
# use dictionary to translate to phyldogID
#  for each speciation node sum number of IES
library(ggtree)
load("~/data/IES_data/rdb/nodeDictionary")
IEScounter <- numeric(7) # count total IES presence probability per speciation node
asr <- read.table("~/data/IES_data/msas/nodeAsr.dat", stringsAsFactors = FALSE, header = TRUE)
geneFamilies <- unique(asr$cluster)
counter <- 1
for(geneFamily in geneFamilies){
  #geneFamily <- geneFamilies[1]
  cat(counter,"/", length(geneFamilies),"\r")
  index <- which(asr$cluster == geneFamily)
  DF <- data.frame(asr[index, ], nodeP = rb2phyldog(asr$node[index], geneFamily), stringsAsFactors = FALSE)
  # find for each phyldog node id what species id it coresponds to
  tr <- read.nhx(paste0("~/data/IES_data/msas/phyldog/results/", geneFamily,".ReconciledTree"))
  dict <- tr@nhx_tags[tr@nhx_tags$Ev == "S", ]
  dict$S <- as.numeric(dict$S)
  # foreach speciation event sum the probability of presence of an IES
  for(S in unique(dict$S)){
#    print(paste(S, sum(DF[DF$nodeP %in% dict$ND[dict$S==S], "presence"])))
    totalP <- sum(DF[DF$nodeP %in% dict$ND[dict$S==S], "presence"])
    if(S == 1 | S == 3 | S == 5 | S == 6){
      if(totalP != trunc(totalP)){
        stop(geneFamily)
      }
    }
    IEScounter[S+1] <-  IEScounter[S+1] + totalP
  }
  counter <- counter + 1
}

knitr::kable(data.frame(spNode = c(0:6), iesNo = IEScounter), digits = 0)
