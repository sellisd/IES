library(seqinr)
path <- "/home/dsellis/data/IES/analysis/msas/aln/"
filteredPath <- "/home/dsellis/data/IES/analysis/msas/filtered/"
singleGenePath <- "/home/dsellis/data/IES/analysis/singleGene/"
if(!dir.exists(filteredPath)){
  dir.create(filteredPath)
}
if(!dir.exists(singleGenePath)){
  dir.create(singleGenePath)
}
files <- dir(path= path,pattern="*.aln.fa$")
nb <- numeric()
id <- numeric()
ps <-  c("PBIA", "PCAU", "POCT", "PPEN", "PPRI", "PSEX", "PSON", "PTET", "PTRE")
pst <- append(ps, "TTHE")
counter <- 1
# single gene gene families (sggf) subset of gene families with one gene from each species (including *T. thermophila*)
sggf <- character()
for(fileName in files){
  print(fileName)
  a <- read.alignment(paste(path,fileName,sep=""),format="fasta")
  nb <- append(nb,a$nb) #number of sequeces
  distM <- dist.alignment(a,matrix="identity")
  id <- append(id,mean(1-distM*distM)) #average pairwise identity
  # check if gene family has one gene from each species
  cur <- sort(substr(a$nam, 0, 4))
  if(identical(cur, ps) || identical(cur, pst)){
    sggf <- append(sggf, fileName)
    counter <- counter + 1
  }  
}
write.table(data.frame(file=files,sequenceNo=nb,avPairId=id),file="~/data/IES/analysis/msas/protAlignStats.dat", quote=FALSE, row.names = FALSE)
highId <- files[which(id>0.5)]
enoughSeq <- files[which(nb>3)]
filtered <- intersect(enoughSeq,highId)
from <- paste(path,filtered,sep="")
to <- paste(filteredPath,filtered,sep="")
r <- file.copy(from, to)
print(all(r)==TRUE)

from <- paste(path, sggf, sep= "")
to <- paste(singleGenePath, sggf, sep = "")
r <- file.copy(from, to)
print(all(r) == TRUE)
