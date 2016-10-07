library(seqinr)
source("~/projects/IES/src/sharedFunctions.R")
source("~/projects/IES/src/sharedPlotFunctions.R")
iess <- read.table("~/data/IES/analysis/mies/silixout/ies.silixout", stringsAsFactors = FALSE, header = FALSE)
d <- read.table("~/data/IES/analysis/mies/miescons.dat", stringsAsFactors = FALSE, header = TRUE)

ppr <- read.table("~/data/IES/analysis/iesdb/ppr.iesdb", header = TRUE, stringsAsFactors = FALSE)
pbi <- read.table("~/data/IES/analysis/iesdb/pbi.iesdb", header = TRUE, stringsAsFactors = FALSE)
pte <- read.table("~/data/IES/analysis/iesdb/pte.iesdb", header = TRUE, stringsAsFactors = FALSE)
ppe <- read.table("~/data/IES/analysis/iesdb/ppe.iesdb", header = TRUE, stringsAsFactors = FALSE)
pse <- read.table("~/data/IES/analysis/iesdb/pse.iesdb", header = TRUE, stringsAsFactors = FALSE)
poc <- read.table("~/data/IES/analysis/iesdb/poc.iesdb", header = TRUE, stringsAsFactors = FALSE)
ptr <- read.table("~/data/IES/analysis/iesdb/ptr.iesdb", header = TRUE, stringsAsFactors = FALSE)
pso <- read.table("~/data/IES/analysis/iesdb/pso.iesdb", header = TRUE, stringsAsFactors = FALSE)
pca <- read.table("~/data/IES/analysis/iesdb/pca.iesdb", header = TRUE, stringsAsFactors = FALSE)
ppr$id <- paste("ppr.", ppr$id, sep = "")
pbi$id <- paste("pbi.", pbi$id, sep = "")
pte$id <- paste("pte.", pte$id, sep = "")
ppe$id <- paste("ppe.", ppe$id, sep = "")
pse$id <- paste("pse.", pse$id, sep = "")
poc$id <- paste("poc.", poc$id, sep = "")
ptr$id <- paste("ptr.", ptr$id, sep = "")
pso$id <- paste("pso.", pso$id, sep = "")
pca$id <- paste("pca.", pca$id, sep = "")
iesdb <- rbind(ppr, pbi, pte, ppe, pse, poc, ptr, pso, pca)
# Construct consensus sequences for candidate mobile IESs

clusters <- c(37918, 35886, 35893, 35880, 402, 12252, 1269, 24810, 27537, 27369, 771, 38004, 27363)
for(clusterId in clusters){
  cluster <- getCoreCluster(clusterId, iess, iesdb)
  write.fasta(sequences = lapply(cluster$seq, s2c), names = cluster$id, file = paste0("~/data/IES/analysis/mies/fasta/", clusterId, ".fa"))
  mafftcmdl <- paste0("mafft --adjustdirectionaccurately --auto /home/dsellis/data/IES/analysis/mies/fasta/", clusterId, ".fa > /home/dsellis/data/IES/analysis/mies/aln/", clusterId,".aln")
  cat(mafftcmdl)
  if (clusterId == clusters[length(clusters)]){
    cat("\n")
  }else{
    cat(" & \n")
  }
}

for(clusterId in clusters){
  weblogocmdl <- paste0("~/tools/weblogo-3.5.0/weblogo -t cl", clusterId, "  -n 200 -f /home/dsellis/data/IES/analysis/mies/aln/", clusterId, ".aln -o /home/dsellis/data/IES/analysis/mies/logo/", clusterId, ".png -F png -A dna")
  cat(weblogocmdl)
  if (clusterId == clusters[length(clusters)]){
    cat("\n")
  }else{
    cat(" & \n")
  }
}

for(clusterId in clusters){
  conscmdl <- paste0("alignbuddy -con /home/dsellis/data/IES/analysis/mies/aln/", clusterId, ".aln |seqbuddy -dm -cs |seqbuddy -ri '.*' cons", clusterId, " > /home/dsellis/data/IES/analysis/mies/fasta/cons", clusterId, ".fa")
  cat(conscmdl)
  if (clusterId == clusters[length(clusters)]){
    cat("\n")
  }else{
    cat(" & \n")
  }
}