# Shared functions
library(seqinr) # make s2c availabel
library(seqLogo)

pairwiseIdentity <- function(a,b){
  # function that returnes a boolean vector with identities of the character matrix with two input strings
  if(length(a) != length(b)){
    stop()
  }
  ac <- toupper(a)
  bc <- toupper(b)
  ac==bc
}

complementarity <- function(a, b){
  if(length(a) != length(b)){
    stop()
  }
  b <- rev(comp(b))
  pairwiseIdentity(a, b)
}

combS <- function(a,b){
  # function that finds all combinations of vectors a and b
  la <- length(a)
  lb <- length(b)
  m <- min(la, lb)
  paired <- matrix(nrow=la*lb,ncol=2)
  counter <- 1
  for(i in c(1:la)){
    for(j in c(1:lb)){
      paired[counter,] <- c(a[i],b[j])
      counter <- counter + 1
    }
  }
  paired
}

prepareLOGO <- function(l){
  aln <- seqinr::as.alignment(nb = length(l),nam=as.character(c(1:length(l))),seq = c(l))
  m <- as.matrix.alignment(aln)
  as <- numeric(0)
  ts <- numeric(0)
  cs <- numeric(0)
  gs <- numeric(0)
  for(i in c(1:length(m[1,]))){
    mt <- table(m[,i])
    as <- append(as,mt["A"])
    cs <- append(cs,mt["C"])
    ts <- append(ts,mt["T"])
    gs <- append(gs,mt["G"])
  }
  #NA to 0
  as[is.na(as)]<-0
  ts[is.na(ts)]<-0
  cs[is.na(cs)]<-0
  gs[is.na(gs)]<-0
  df <- data.frame("A"=as,"C"=cs,"G"=gs,"T"=ts)
  seqLogo(makePWM(t(df/rowSums(df))))
}

#remove NAs and names from vectors and
nadel <- function(x){
  as.vector(x[!is.na(x)])
}
