library(seqinr)
path <- "/home/dsellis/data/IES_data/msas/alignments/aln/"
filteredPath <- "/home/dsellis/data/IES_data/msas/alignments/filtered/"
if(!dir.exists(filteredPath)){
  dir.create(filteredPath)
}
files <- dir(path= path,pattern="*.aln$")
nb <- numeric()
id <- numeric()
for(file in files){
  print(file)
  a <- read.alignment(paste(path,file,sep=""),format="clustal")
  nb <- append(nb,a$nb) #number of sequeces
  id <- append(id,mean(dist.alignment(a,matrix="identity"))) #average pairwise identity
}
write.table(data.frame(file=files,sequenceNo=nb,avPairId=id),file="~/data//IES_data/msas/protAlignStats.dat", quote=FALSE)
highId <- files[which(id>0.5)]
enoughSeq <- files[which(nb>3)]
filtered <- intersect(enoughSeq,highId)
from <- paste(path,filtered,sep="")
to <- paste(filteredPath,filtered,sep="")
r <- file.copy(from, to)
print(all(r)==TRUE)
