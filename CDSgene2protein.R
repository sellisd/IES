# read cds coordinates in gene coordinates and calculate the corresponding protein coordinates

ptecds <- read.table("~/data/IES_data/ptetraurelia/pte.cds", header = TRUE, stringsAsFactors = FALSE)
pbicds <- read.table("~/data/IES_data/pbiaurelia/pbi.cds", header = TRUE, stringsAsFactors = FALSE)
ptecdsprot <- gene2protCDS(ptecds)
pbicdsprot <- gene2protCD(pbicds)
write.table(ptecdsprot, file = "~/data/IES_data/ptetraurelia/pte.cds.prot", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(pbicdsprot, file = "~/data/IES_data/pbiaurelia/pbi.cds.prot", row.names = FALSE, quote = FALSE, sep = "\t")
