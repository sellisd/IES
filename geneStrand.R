# find strand of genes
# --- P. tetraurelia ---
pteG <- read.table("~/data/IES_data/ptetraurelia/ptetraurelia_mac_51_annotation_v2.0.gff3", as.is = TRUE)
geneIds <- pteG[pteG$V3 == "gene", 9]
matched <- regexpr("ID=[^;]*;", geneIds, perl = TRUE)
# clean names fron ID=...;
geneIds <- gsub("(ID=)|;", "", regmatches(geneIds, matched), perl = TRUE)
# G to P
geneIds <- sub("G","P",geneIds)
pteStrand <- pteG$V7[pteG$V3 == "gene"]
names(pteStrand) <- geneIds
# --- P. biaurelia ---
pbiG <- read.table("~/data/IES_data/pbiaurelia/pbiaurelia_V1-4_annotation_v2.0.gff3", as.is = TRUE)
geneIds <- pbiG[pbiG$V3 == "gene", 9]
matched <- regexpr("ID=[^;]*;", geneIds, perl = TRUE)
geneIds <- gsub("(ID=)|;", "", regmatches(geneIds, matched), perl = TRUE)
geneIds <- sub("G","P", geneIds)
pbiStrand <- pbiG$V7[pbiG$V3 == "gene"]
names(pbiStrand) <- geneIds

save(pteStrand, pbiStrand, file = "~/data/IES_data/rdb/geneStrand")
