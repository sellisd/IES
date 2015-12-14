#calculate transitions from nodes to nodes
library(tidyr)
library(reshape2)
source("~/projects/IES/src/sharedFunctions.R")
burnIn <- 1000
load("~/data/IES_data/rdb/charMats")
load("~/data/IES_data/rdb/nodeDictionary")
# for all gene families with IES in conserved blocks
ancestralStates <- data.frame(cluster = numeric(0),
                              Iteration = numeric(0),
                              rb = numeric(0),
                              iesColumn = numeric(0),
                              presenceAbsence = numeric(0),
                              r = numeric(0))

counter <- 1
clusters <- unique(charMats$cluster)
for(cluster in clusters){
  #cluster <- clusters[1]
  print(paste(counter, "/", length(clusters)))
  # parse ancestral rate reconstructions to a R-friendly format
  # read ancestral state reconstructions and transform it to long data.frame
  ancStates <- read.table(paste("~/data/IES_data/msas/asr/run4/ancStates", cluster, ".log", sep = ""), header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  tidyAncSt <- gather(ancStates, rb, paColumns, -Iteration)
  # how many columns should we expect? check the first row and treat it as a character (in case there is only one)
  iesColumnNumber <- length(strsplit(as.character(tidyAncSt$paColumns[1]),",", fixed = TRUE)[[1]])
  tidyAncSt <- separate(tidyAncSt, paColumns, into = c(1:iesColumnNumber), sep = ",", fixed = TRUE)
  tidyAncSt <- melt(tidyAncSt, id.vars = c("Iteration", "rb"), variable.name = "iesColumn", value.name = "presenceAbsence")
  # remove burnIn
  tidyAncSt <- tidyAncSt[tidyAncSt$Iteration > burnIn, ]
  # make numeric presence Absence so we can calculate means and stdevs
  tidyAncSt$presenceAbsence <- as.numeric(tidyAncSt$presenceAbsence)
  tidyAncSt$rb <- substr(tidyAncSt$rb, 5, length(tidyAncSt$rb)) # rename nodes from end_0 to 0
  # translate numbering to R format
  tidyAncSt <- data.frame(cluster = cluster, tidyAncSt, r = rb2r(tidyAncSt$rb, cluster = cluster), stringsAsFactors = FALSE)
  counter <- counter + 1
}
save(tidyAncSt, file = "~/data/IES_data/rdb/ancestralStates")

