# find IES that were recently gained
library(seqinr)  
library(dplyr)

load("~/data/IES_data/rdb/charMats")

source("~/projects/fgmo/colors.R")
source("~/projects/IES/src/sharedFunctions.R")
source("~/projects/IES/src/sharedPlotFunctions.R")

pdf("~/projects/IES/reports/26.IESLoss/recentGLfiltered.pdf")
par(mfcol=c(2,2), las = 1)
# filter out small scaffolds
# load scaffold lengths
pbiSc <- read.table("~/data/IES/analysis/pbi.scaf", stringsAsFactors = FALSE, header = TRUE)
pteSc <- read.table("~/data/IES/analysis/pte.scaf", stringsAsFactors = FALSE, header = TRUE)
largeScafs <- c(pbiSc$scaffold[pbiSc$length>10000], pteSc$scaffold[pteSc$length>10000])

#filter only events from speciation node 4 to 5 or 4 to 6
gl <- read.table("~/data/IES_data/msas/recentGainLossIES.dat", header = TRUE, stringsAsFactors = FALSE)

recentGl <- gl[which(gl$fromS == "4" & (gl$toS == "5" | gl$toS == "6")),]

# keep only events that are in large scaffolds
pref <- nrow(recentGl) #number of events before filtering
recentGl <- recentGl[geneInScaffold(recentGl$gene) %in% largeScafs, ]
postf <- nrow(recentGl)

print(paste(pref - postf, "IES in small scaffolds"))

# filter out events close to each other, close to CDS boundaries and close to Gblocks boundaries
clusters <- unique(gl$cluster)
eventCounter <- data.frame(pbiLoss = 0, pbiGain = 0, pteLoss = 0, pteGain = 0, stringsAsFactors = FALSE)
clusterCounter <- 1
total <- length(clusters) 
for(cluster in clusters){
  print(paste0(clusterCounter,"/", total))
  clusterCounter <- clusterCounter + 1
  eventsInClusterI <- which(recentGl$cluster == cluster)
  events <- recentGl[eventsInClusterI, c("iesColumn", "eventType", "gene")]
  #keep only events from the ancestor of *P. tetraurelia* and *P. biaurelia* to either *P. tetraurelia* or *P. biaurelia* (from speciation node 4 to speciation nodes 5 or 6)
  if(nrow(events) == 0){ # if we have no events
    next
  }
  events$iesColumn <- events$iesColumn + 1 # from 0-indexed to 1-indexed
  # find gblock boundaries
  gblocks <- read.table(paste("~/data/IES_data/msas/alignments/filtered/cluster.", cluster, ".aln.fasta.gblocks", sep = ""), stringsAsFactors = FALSE)
  # find CDS bounaries
  cdsboundAll <- read.table(paste("~/data/IES_data/msas/alignments/filtered/cluster.", cluster, ".cds.boundaries", sep = ""), stringsAsFactors = FALSE)
  pad <- 4
  gblocksCenterIR <- IRanges(start = gblocks$V2 + pad, end = gblocks$V3 - pad)
  #if CDS is too small remove it completely
  cdsbound <- cdsboundAll[which(cdsboundAll$V2+pad < cdsboundAll$V3-pad), ]
  cdsCenterIR <- IRanges(start = cdsbound$V2 + pad, end = cdsbound$V3 - pad)
  loc <- unique(charMats[which(charMats$cluster == cluster), c("column", "begin", "end")])
  locIR <- IRanges(start = loc$begin, end = loc$end)
  extEvents <- left_join(events, loc, by = c("iesColumn" = "column"))
  extEventsIR <- IRanges(start = extEvents$begin, end = extEvents$end)
  # keep events that are not close to the borders + padding
  inblocks <- overlapsAny(extEventsIR, gblocksCenterIR, type = "within")
  incds <- overlapsAny(extEventsIR, cdsCenterIR, type = "within")
  keepIndex <- inblocks & incds
  filtered <- extEvents[keepIndex, ]
  #filtered <- extEvents[inblocks, ]
  if(nrow(filtered) == 0){
    next # no events that are not close to a block boundary
  }
  # check if we have both gains and losses in the same alignment
  if(length(unique(filtered$eventType)) == 2){
    # remove cases where gains and losses are too close to each other
    gains <- filtered[filtered$eventType == "gain", ]
    losses <- filtered[filtered$eventType == "loss", ]
    gainsIR <- IRanges(start = gains$begin, end = gains$end)
    lossesIR <- IRanges(start = losses$begin, end = losses$end)
    d2n <- distanceToNearest(gainsIR, subject = lossesIR)
    if(min(d2n@elementMetadata$distance) < pad){
      next
      # ignore gene families that have both gains and losses
    }
  }
  sp <- gene2species(filtered$gene)
  eventCounter$pbiLoss <- eventCounter$pbiLoss + sum(filtered$eventType == "loss" & sp == "Paramecium_biaurelia")
  eventCounter$pteLoss <- eventCounter$pteLoss + sum(filtered$eventType == "loss" & sp == "Paramecium_tetraurelia")
  eventCounter$pbiGain <- eventCounter$pbiGain + sum(filtered$eventType == "gain" & sp == "Paramecium_biaurelia")
  eventCounter$pteGain <- eventCounter$pteGain + sum(filtered$eventType == "gain" & sp == "Paramecium_tetraurelia")
  alignF <- paste0("/home/dsellis/data/IES_data/msas/alignments/filtered/cluster.", cluster, ".nucl.fa")
  aln <- read.alignment(file = alignF, format = "fasta")
  alnM <- as.matrix.alignment(aln)
  introns <- data.frame(gene = gene2protName(cdsboundAll$V1), begin = cdsboundAll$V2, end = cdsboundAll$V3, stringsAsFactors = FALSE)
  for(i in c(1:nrow(filtered))){
    printIESalign(geneFamily = cluster, iesColumn = filtered$iesColumn[i], pad = 5, filtered[i, ], introns)
    title(paste(cluster, filtered$iesColumn[i], filtered$eventType[i]))
  }
}
dev.off()