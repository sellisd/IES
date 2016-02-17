# extract intron information from gff3 files to calculate some descriptive statistics
library(parallel)
source("~/projects/IES/src/sharedFunctions.R")
ppr <- read.table("~/data/IES/analysis/ppr.cds", stringsAsFactors = FALSE, header = TRUE)
pbi <- read.table("~/data/IES/analysis/pbi.cds", stringsAsFactors = FALSE, header = TRUE)
pte <- read.table("~/data/IES/analysis/pte.cds", stringsAsFactors = FALSE, header = TRUE)
ppe <- read.table("~/data/IES/analysis/ppe.cds", stringsAsFactors = FALSE, header = TRUE)
pse <- read.table("~/data/IES/analysis/pse.cds", stringsAsFactors = FALSE, header = TRUE)
poc <- read.table("~/data/IES/analysis/poc.cds", stringsAsFactors = FALSE, header = TRUE)
ptr <- read.table("~/data/IES/analysis/ptr.cds", stringsAsFactors = FALSE, header = TRUE)
pca <- read.table("~/data/IES/analysis/pca.cds", stringsAsFactors = FALSE, header = TRUE)

cds <- list(ppr, pbi, pte, ppe, pse, poc, ptr, pca)
intronStats <- mclapply(cds, introns, mc.cores = length(cds), mc.silent = TRUE)
rdb <- "~/data/IES/analysis/rdb"
if(!dir.exists(rdb)){
  dir.create(rdb)
}
save(intronStats, file = file.path(rdb,"intronStats"))
