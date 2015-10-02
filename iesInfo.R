# find IES info

# find floating IES and IES border sequence
windowSize <- 20
# lengthClass
# ------- for P. tetraurelia----------
ptemacseq <- read.table("~/data/IES_data/ptetraurelia/Pte.ies.mac_seq", as.is=TRUE, row.names = 1)
upstream <- substr(ptemacseq$V2,0,15)
downstream <- substr(ptemacseq$V2,18,32)
flanks <- data.frame(upstream = upstream, downstream = downstream, row.names = row.names(ptemacseq), stringsAsFactors = FALSE)

pte <- read.table("~/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3", as.is = TRUE)
pteL <- strsplit(pte[,9],";",fixed=TRUE)
iesNo <- nrow(pte)
id <- character(iesNo)
front <- character(iesNo)
back <- character(iesNo)
IESlength <- numeric(iesNo)
IESsequence <- character(iesNo)
start <- pte$V4
end <- pte$V5
upstream <- character(iesNo)
downstream <- character(iesNo)
floatingB <- (end-start>1)
counter <- 1
for(i in pteL){
  sbL <- strsplit(i,"=")
  for(j in sbL){
    if(j[1]=="ID"){
      id[counter] <- j[2]
      upstream[counter] <- flanks[j[2],"upstream"]
      downstream[counter] <- flanks[j[2],"downstream"]
    }
    if(j[1]=="sequence"){
      IESseq <- j[2]
      l <- nchar(IESseq)
      front[counter] <- substr(IESseq,0,windowSize)
      back[counter] <- substr(IESseq,l-windowSize+1,l)
      IESsequence[counter] <- IESseq
      IESlength[counter] <- l-2 # The IES length does not include both TAs
    }
  }
  counter <- counter + 1
}
pteD <- data.frame(length=IESlength, seq = IESsequence, front=front,back=back, start = start, end = end, floating = floatingB, upstream = upstream, downstream = downstream, row.names=id, stringsAsFactors=FALSE)

# ------- for P. biaurelia ----------
pbi <- read.table("~/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.fl.gff3", as.is = TRUE)
pbiL <- strsplit(pbi[,9],";",fixed=TRUE)
iesNo <- nrow(pbi)
id <- character(iesNo)
front <- character(iesNo)
back <- character(iesNo)
IESlength <- numeric(iesNo)
IESsequence <- character(iesNo)
start <- pbi$V4
end <- pbi$V5
floatingB <- (end-start>1)
upstream <- character(iesNo)
downstream <- character(iesNo)
counter <- 1
for(i in pbiL){
  sbL <- strsplit(i,"=")
  for(j in sbL){
    if(j[1] == "ID"){
      id[counter] <- j[2]
    }else if(j[1] == "sequence"){
      IESseq <- j[2]
      l <- nchar(IESseq)
      front[counter] <- substr(IESseq,0,windowSize)
      back[counter] <- substr(IESseq,l-windowSize+1,l)
      IESsequence[counter] <- IESseq
      IESlength[counter] <- l-2 # The IES length does not include both TAs
    }else if(j[1] == "mac_seq"){
      upstream[counter] <- substr(j[2], 0, 15)
      downstream[counter] <- substr(j[2],18,32)
    }
  }
  counter <- counter + 1
}
pbiD <- data.frame(length = IESlength, seq = IESsequence, front = front, back = back, start = start, end = end, floating = floatingB, upstream = upstream, downstream = downstream, row.names=id, stringsAsFactors=FALSE)
save(pteD, pbiD, file = "~/data/IES_data/rdb/iesInfo")
