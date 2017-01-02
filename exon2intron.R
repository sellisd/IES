# read exon bed files and print corresponding intron bed files
setwd('.')
source("sharedFunctions.R")
args<-commandArgs(TRUE);
inFile <- args[1]
outFile <- args[2]
#Rscript --vanilla exon2intron.R inFile outFile
exons <- read.table(inFile, stringsAsFactors = FALSE)
names(exons) <- c("scaffold", "start", "end", "name", "gene")
introns <- exon2intronsBed(exons)
write.table(introns, file = outFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)