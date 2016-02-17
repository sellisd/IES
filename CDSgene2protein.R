# read cds coordinates in gene coordinates and calculate the corresponding protein coordinates
source("~/projects/IES/src/sharedFunctions.R")
ptecds <- read.table("~/data/IES_data/ptetraurelia/pte.cds", header = TRUE, stringsAsFactors = FALSE)
ptecdsprot <- gene2protCDS(ptecds)
write.table(ptecdsprot, file = "~/data/IES_data/ptetraurelia/pte.cds.prot", row.names = FALSE, quote = FALSE, sep = "\t")

pbicds <- read.table("~/data/IES_data/pbiaurelia/pbi.cds", header = TRUE, stringsAsFactors = FALSE)
pbicdsprot <- gene2protCDS(pbicds)
write.table(pbicdsprot, file = "~/data/IES_data/pbiaurelia/pbi.cds.prot", row.names = FALSE, quote = FALSE, sep = "\t")
