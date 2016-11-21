#
# filter BLAST HSPs
#
blout <- read.table("/home/dsellis/data/IES/analysis/mies/blastout/ies.blastout", stringsAsFactors = FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

source("~/projects/IES/src/sharedFunctions.R")
#find length of all IESs from qseqid and apply goodHSP
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


iesl <- iesdb$length[match(blout$qseqid, iesdb$id)]
iesl <- iesl + 2 # add 2 for the extra TA
# each IES should have 1*** AND *1** OR **1* AND ***1
f <- logical(nrow(blout))
iesIds <- unique(blout$qseqid)
counter <- 1
for(iesId in iesIds){
  cat(counter, "/", length(iesIds), "\r")
  I <- which(blout$qseqid == iesId)
  hspTypes <- HSPtype(rep(20, length(I)), blout$qstart[I], blout$qend[I], iesl[I])
  if((any(hspTypes[,1] == TRUE) & any(hspTypes[,2] == TRUE)) | (any(hspTypes[,3] == TRUE) & any(hspTypes[,4] == TRUE))){
    f[I] <- TRUE
  }
  counter <- counter + 1
}

# and compare BLAST results

nonselfHits <- blout$qseqid != blout$sseqid
iesbl <- blout[f&nonselfHits,]
iesflblout <- read.table("/home/dsellis/data/IES/analysis/mies/blastout/iesflanks.blastout", stringsAsFactors = FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
qsi <- paste0(iesbl$qseqid, iesbl$sseqid)
qsf <- paste0(iesflblout$qseqid, iesflblout$sseqid)
sqf <- paste0(iesflblout$sseqid, iesflblout$qseqid)

hl <- (qsi %in% qsf) | (qsi %in%sqf)
# homologous loci (conserved IESs)
write.table(iesbl[hl, ], file = "~/data/IES/analysis/mies/blastout/ies.homl.blastout", quote = FALSE, col.names = FALSE, row.names = FALSE)
# non-homologous loci (mobile IESs)
write.table(iesbl[!hl, ], file = "~/data/IES/analysis/mies/blastout/ies.nonhoml.blastout", quote = FALSE, col.names = FALSE, row.names = FALSE)
