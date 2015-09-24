library(seqinr)
path <- "/home/dsellis/data/IES_data/msas/alignments/aln/"
filteredPath <- "/home/dsellis/data/IES_data/msas/alignments/filtered/"
if(!dir.exists(filteredPath)){
  dir.create(filteredPath)
}
files <- dir(path= path,pattern="*.aln$")
nb <- numeric()
id <- numeric()
for(fileName in files){
  print(fileName)
  a <- read.alignment(paste(path,fileName,sep=""),format="clustal")
  nb <- append(nb,a$nb) #number of sequeces
  distM <- dist.alignment(a,matrix="identity")
  id <- append(id,mean(1-distM*distM)) #average pairwise identity
}
write.table(data.frame(file=files,sequenceNo=nb,avPairId=id),file="~/data/IES_data/msas/protAlignStats.dat", quote=FALSE)
highId <- files[which(id>0.5)]
enoughSeq <- files[which(nb>3)]
filtered <- intersect(enoughSeq,highId)
from <- paste(path,filtered,sep="")
to <- paste(filteredPath,filtered,sep="")
r <- file.copy(from, to)
print(all(r)==TRUE)
