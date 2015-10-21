# for each species
source("~/projects/IES/src/sharedFunctions.R")
pbiG <- read.table("~/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3", as.is = TRUE)
pteG <- read.table("~/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3", as.is = TRUE)
pseG <- read.table("~/data/IES_data/psexaurelia/psexaurelia_AZ8-4_annotation_v2.0.gff3", as.is = TRUE)
pcaG <- read.table("~/data/IES_data/pcaudatum_43c3d_annotation_v2.0/pcaudatum_43c3d_annotation_v2.0.gff3", as.is = TRUE)

# extract information for genes and order data.frames
cat("extracting gene information\n")
pbiBE <- extractGenes(pbiG)
pteBE <- extractGenes(pteG)
pseBE <- extractGenes(pseG)
pcaBE <- extractGenes(pcaG)

# find closest neighbor of each gene
cat("finding neighbors in P. biaurelia\n")
pbiGP <- genePairs(pbiBE)
cat("finding neighbors in P. tetraurelia\n")
pteGP <- genePairs(pteBE)
cat("finding neighbors in P. sexaurelia\n")
pseGP <- genePairs(pseBE)
cat("finding neighbors in P. caudatum\n")
pcaGP <- genePairs(pcaBE)
cat("saving results\n")
save(pbiGP, pteGP, pseGP, pcaGP, file = "~/data/IES_data/intergenic/genePairs")
cat("done!\n")