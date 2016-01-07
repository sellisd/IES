load("~/data/IES_data/rdb/nodePairs")
source("~/projects/IES/src/sharedFunctions.R")
#from pairs of nodes and ancestral states DF count transitions
gainLoss <- data.frame(cluster = numeric(0),
                        iesColumn = numeric(0),
                        fromEvent = numeric(0),
                        toEvent = numeric(0),
                        fromR = numeric(0),
                        toR = numeric(0),
                        presence = numeric(0),
                        absence = numeric(0),
                        gain = numeric(0),
                        loss = numeric(0))
clusters <- read.table("~/data/IES_data/msas/asr/geneFamilies.dat", stringsAsFactors = FALSE, as.is = TRUE)
gfCounter <- 1
for(cluster in clusters$V1[1:1000]){
  print(paste(gfCounter, "/", length(clusters$V1)))
  rFileF <- paste0("~/data/IES_data/rdb/ancestralStates/", cluster)
  if(!file.exists(rFileF)){
    stop(paste(rFileF), " is missing")
  }
  load(rFileF) 
  clusterI <- which(nodePairs$cluster == cluster)
  for(iesColumn in unique(ancestralStates$iesColumn)){
    #iesColumn <- 1
    columnIndex <- which(ancestralStates$iesColumn == iesColumn)
    # for each node Pair calculate events
    l <- nrow(nodePairs[clusterI, ])
    presence <- numeric(l)
    absence <- numeric(l)
    gain <- numeric(l)
    loss <- numeric(l)
    counter <- 1
    for(i in clusterI){
      from <-  nodePairs[i,"fromR"]
      to <- nodePairs[i,"toR"]
      fromIndex <- which(ancestralStates$r == from)
      toIndex <- which(ancestralStates$r == to)
      fromDF <- ancestralStates[intersect(columnIndex, fromIndex),]
      toDF <- ancestralStates[intersect(columnIndex, toIndex),]
      fromToDF <- merge(fromDF, toDF, by = "Iteration")
      et <- eventType(fromToDF$presenceAbsence.x, fromToDF$presenceAbsence.y)
      presence[counter] <- et[1]
      absence[counter] <- et[2]
      gain[counter] <- et[3]
      loss[counter] <- et[4]
      counter <- counter + 1
    }
    gainLoss <- rbind(gainLoss, data.frame(cluster = cluster,
                                           iesColumn = iesColumn,
                                           fromEvent = nodePairs$fromEvent[clusterI],
                                           toEvent = nodePairs$toEvent[clusterI],
                                           fromR = nodePairs$fromR[clusterI],
                                           toR = nodePairs$toR[clusterI],
                                           presence = presence,
                                           absence = absence,
                                           gain = gain,
                                           loss = loss,
                                           stringsAsFactors = FALSE)
    )
  }
  gfCounter <- gfCounter + 1
}
save(gainLoss, file = "~/data/IES_data/rdb/gainLoss")
