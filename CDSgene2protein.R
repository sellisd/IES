# read cds coordinates in gene coordinates and calculate the corresponding protein coordinates
args<-commandArgs(TRUE);
cdsFile <- args[1]
outFile <- args[2]
#cdsFile "~/data/IES_data/ptetraurelia/pte.cds"
#outFile "~/data/IES_data/ptetraurelia/pte.cds.prot"
#Rscript --vanilla CDSgene2protein.R cdsFile outFile
source("~/projects/IES/src/sharedFunctions.R")
cds <- read.table(cdsFile, header = TRUE, stringsAsFactors = FALSE)
cdsprot <- gene2protCDS(cds)
write.table(cdsprot, file = outFile, row.names = FALSE, quote = FALSE, sep = "\t")
