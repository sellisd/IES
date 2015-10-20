# Shared functions

# required libraries
library(seqinr) # make s2c availabel
library(seqLogo)
library(plyr)

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
  # calculate per locus complementarity of two sequences, by reverse complementing the second and calculating identity
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
# 
# x <- data.frame(letters = letters[1:10], numbers1 = c(1:10), numbers2 = c(1,1,1,1,1,1,1,2,2,2), ignored = letters[1:10], numbers3 = c(NA,NA,NA,2,3,3,4,4,5,5))
# col2permute <- c("letters", "numbers1")
# col2keepFixed <- c("numbers2", "numbers3")

permuteBy <- function(x, col2permute, col2keepFixed, ret = "index"){
  # function that permutes values in data.frame column but keeping values in another column fixed
  # col2permute can be a column name of the x data.frame  or a character vector with multiple such names
  # col2keepFixed can be a column name of the x data.frame
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
  permutedI <- numeric(nrow(x))
  #make a key from all columns that we need to keep
  if(length(col2keepFixed) >  1){
    key <-  apply(x[ , col2keepFixed], 1, paste0, collapse =" ")
  }else{
    key <- x[, col2keepFixed]
  }
  #find unique combinations
  lev <- unique(key)
  #find index of unique combinations
  for(l in lev){
    index <- which(key == l)
    if(length(index) == 1){
      newIndex = index
    }else{
      newIndex <- sample(index, length(index))    
    }
    permutedI[index] <- newIndex
  }
  for(i in col2permute){
    x[i] <- x[permutedI, i]
  }
  if(ret == "index"){
    return(permutedI)
  }else if(ret == "data.frame"){
    return(x)
  }else{
    stop("return error")
  }
}

compl <- function(gc){
  # The expected complementarity of random sequences with gc GC%
  at <- 1 - gc  
  at*at/2+gc*gc/2
}

appa <- function(x){
  # approximate presence-absence
  # The function determines the state of an IES in a species from character vector containing information from all paralogs
  # INPUT: a character string with the following options
  #       0:  absence
  #       NA: unknown IES state
  #       !=0: presence
  # OUTPUT: NA: no IES information
  #          0: all absent
  #          1: at least one present
  #         -1: empty input
  if(length(x) == 0){
    return(-1)
  }
  if(all(is.na(x))){
    return(NA)
  }
  if(all(x == 0, na.rm = TRUE)){
    return(0)
  }
  if(any(x != 0)){
    return(1)
  }  
}
#sp <- spGT
presenceAbsence <- function(sp){
  # Function to assign patterns of presence / absence of IES in subtrees.
  # In the case of multiple paralogs, presence is counted as at least one parlaog having an IES
  load("~/data/IES_data/rdb/charMats")
  # add to charMats a column with in-/outgroup/NA and species
  extMat <- join(charMats, sp, by = "geneId")
  biI <- which(extMat[, "species"] == "Paramecium_biaurelia")
  teI <- which(extMat[, "species"] == "Paramecium_tetraurelia")
  seI <- which(extMat[, "species"] == "Paramecium_sexaurelia")
  caI <- which(extMat[, "species"] == "Paramecium_caudatum")
  uniqueSubtrees <- unique(extMat$subtree)
  # the total number of patterns will be the product of the number of subtrees by the number of columns in each one
  countRows <- extMat[,c("column","subtree")]
  countRows <- countRows[!is.na(countRows$subtree), ]
  countRows <- unique(paste(countRows[,1], countRows[,2]))
  patternsM <- matrix(ncol = 4, nrow = length(countRows)) # minimum length of matrix
#  patternsV <- character(length(uniqueSubtrees))
  counter <- 1
  #for(clustcol in unique(homIES)){
  for(subtree in uniqueSubtrees){
    #subtree <- uniqueSubtrees[2]
    if(is.na(subtree)){
      next
    }
    indexCM <- which(extMat$subtree == subtree) # create and index for the current subtree
    homologIES <- unique(extMat[indexCM, "column"])
    for(column in homologIES){
      # column <- homologIES[1]                 # and one for each column
      indexCMC <- which(extMat$subtree == subtree & extMat$column == column)
      # for each species find if at least one IES
      patternsM[counter, ] <- c(appa(extMat[intersect(biI, indexCMC), "ies"]),
                                appa(extMat[intersect(teI, indexCMC), "ies"]),
                                appa(extMat[intersect(seI, indexCMC), "ies"]),
                                appa(extMat[intersect(caI, indexCMC), "ies"]))
      counter <- counter + 1
    }
  }
  patternsM
}


gene2species <- function(string){
  # get a character vector of species names from 
  # a character vector of gene names
  knownSpeciesNames <- c("Paramecium_caudatum", "Paramecium_sexaurelia", "Paramecium_tetraurelia", "Paramecium_biaurelia", "Tetrahymena_thermophila")
  names(knownSpeciesNames) <-  c("PCAU", "PSEX", "PTET", "PBIA", "TTHE")
  speciesNames <- knownSpeciesNames[substr(string, 0, 4)]
  if(any(is.na(speciesNames))){
    stop(paste("unknown name", string[which(is.na(speciesNames))]))
  }
  speciesNames
}

extractGenes <- function(ggff){
  # extract information for genes from a gff3 file
  geneIndex <- which(ggff$V3 == "gene")
  geneIds <- ggff[geneIndex, 9]
  matched <- regexpr("ID=[^;]*;", geneIds, perl = TRUE)
  # clean names fron ID=...;
  geneIds <- gsub("(ID=)|;", "", regmatches(geneIds, matched), perl = TRUE)
  # G to P 
  geneIds <- sub("G","P",geneIds)
  scaffold <- ggff[geneIndex, 1]
  strand <- ggff[geneIndex, 7]
  begin <- ggff[geneIndex, 4]
  end <- ggff[geneIndex, 5]
  df <- data.frame(id = geneIds, begin = begin, end = end, strand = strand, scaffold = scaffold, stringsAsFactors = FALSE)
  #order by scaffold, then by start and then by end
  df[order(df$scaffold, df$begin, df$end), ]
}


overlapping <- function(a, b){
  # function that tests if two elements a and b are overlapping
  if(a$strand == '-'){
    beginA <- a$end
    endA <- a$begin
  }else{
    beginA <- a$begin
    endA <- a$end
  }
  if(b$strand == '-'){
    beginB <- b$end
    endB <- b$begin
  }else{
    beginB <- b$begin
    endB <- b$end
  }
  if(beginB >= beginA & beginB <= endA){  # B starts within A
    return(1) 
  }
  if(endB >= beginA & endB <= endA){      # B ends within A
    return(1)
  }
  if(endB >= endA & beginB <= beginA){
    return(1)                             # A is within B
  }
  return(0)
  # # testing code for function
  # a <- data.frame(begin = 10, end = 20, strand = "+") 
  # b1 <- data.frame(begin = 1, end = 2, strand = "+")
  # b2 <- data.frame(begin = 1, end = 10, strand = "+")
  # b3 <- data.frame(begin = 1, end = 12, strand = "+")
  # b4 <- data.frame(begin = 1, end = 20, strand = "+")
  # b5 <- data.frame(begin = 1, end = 22, strand = "+")
  # b6 <- data.frame(begin = 10, end = 20, strand = "+")
  # b7 <- data.frame(begin = 12, end = 18, strand = "+")
  # b8 <- data.frame(begin = 12, end = 20, strand = "+")
  # b9 <- data.frame(begin = 12, end = 25, strand = "+")
  # b10 <- data.frame(begin = 20, end = 25, strand = "+")
  # b11 <- data.frame(begin = 22, end = 25, strand = "+")
  # 
  # overlapping(a, b1) # 0
  # overlapping(a, b2) # 1
  # overlapping(a, b3) # 1
  # overlapping(a, b4) # 1
  # overlapping(a, b5) # 1
  # overlapping(a, b6) # 1
  # overlapping(a, b7) # 1
  # overlapping(a, b8) # 1
  # overlapping(a, b9) # 1
  # overlapping(a, b10) # 1
  # overlapping(a, b11) # 0
  
}
