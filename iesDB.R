#idsDB
# make a database of commonly used IES information for all species in a handy (and uniform) format
#
# INPUT: 
#   ~/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.gff3
#   ~/data/IES_data/pbiaurelia//internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.ies
#   ~/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.ies
#   ~/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3
#
# OUTPUT: 
#  R saved data.frames, one for each species with: length front back and id as row.names
#
# parameters: 
# windowSize <- 10, window size of borders of IES sequence retained
windowSize <- 10
#-----P. caudatum-----------
pca <- read.table("~/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.gff3",as.is=TRUE)
pcaL <- strsplit(pca[,9],";",fixed=TRUE)
id <- character()
front <- character()
back <- character()
IESlength <- numeric()

for(i in pcaL){
  sbL <- strsplit(i,"=")
  for(j in sbL){
    if(j[1]=="ID"){
      id <- append(id,j[2])
    }
    if(j[1]=="sequence"){
      IESseq <- j[2]
      l <- nchar(IESseq)
      IESlength <- append(IESlength,l)
      front <- append(front,substr(IESseq,0,windowSize))
      # P caudatum does not have the ending TA
      back <- append(back,paste(substr(IESseq,l-windowSize+1+2,l),"TA",sep=""))
    }
  }
}
pcaD <- data.frame(length=IESlength,front=front,back=back,row.names=id, stringsAsFactors=FALSE)

#--- P. biaurelia, P. sexuarelia-------
pbi <- read.table("~/data/IES_data/pbiaurelia//internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.ies", as.is = TRUE)
pbiL <- strsplit(pbi[,6],";",fixed=TRUE)

IESseq <- pbi$V5
IESlength  <- pbi$V3
front <- substr(IESseq,0,windowSize)
back <- substr(IESseq,IESlength+2-windowSize+1,IESlength+2)
id <- character()
for(i in pbiL){
  if(substr(i[1],0,3)=="ID="){
    id <- append(id,substr(i[1],4,nchar(i[1])))
  }
}
pbiD <- data.frame(length=IESlength,front=front,back=back,row.names=id, stringsAsFactors=FALSE)

pse <- read.table("~/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.ies", as.is = TRUE)
pseL <- strsplit(pse[,6],";",fixed=TRUE)

IESseq <- pse$V5
IESlength  <- pse$V3
front <- substr(IESseq,0,windowSize)
back <- substr(IESseq,IESlength+2-windowSize+1,IESlength+2)
id <- character()
for(i in pseL){
  if(substr(i[1],0,3)=="ID="){
    id <- append(id,substr(i[1],4,nchar(i[1])))
  }
}
pseD <- data.frame(length=IESlength,front=front,back=back,row.names=id, stringsAsFactors=FALSE)

#------P. tetraurelia-----------
pte <- read.table("~/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3", as.is = TRUE)
pteL <- strsplit(pte[,9],";",fixed=TRUE)
id <- character()
front <- character()
back <- character()
IESlength <- numeric()

for(i in pteL){
  sbL <- strsplit(i,"=")
  for(j in sbL){
    if(j[1]=="ID"){
      id <- append(id,j[2])
    }
    if(j[1]=="sequence"){
      IESseq <- j[2]
      l <- nchar(IESseq)
      front <- append(front,substr(IESseq,0,windowSize))
      back <- append(back,substr(IESseq,l-windowSize+1,l))
      IESlength <- append(IESlength,l-2) # The IES length does not include both TAs
    }
  }
}
pteD <- data.frame(length=IESlength,front=front,back=back,row.names=id, stringsAsFactors=FALSE)

save(pcaD,pbiD,pteD,pseD,file="~/data/IES_data/iesDB")
