
# test function overlapping

a <- data.frame(begin = 10, end = 20, strand = "+") 
b1 <- data.frame(begin = 1, end = 2, strand = "+")
b2 <- data.frame(begin = 1, end = 10, strand = "+")
b3 <- data.frame(begin = 1, end = 12, strand = "+")
b4 <- data.frame(begin = 1, end = 20, strand = "+")
b5 <- data.frame(begin = 1, end = 22, strand = "+")
b6 <- data.frame(begin = 10, end = 20, strand = "+")
b7 <- data.frame(begin = 12, end = 18, strand = "+")
b8 <- data.frame(begin = 12, end = 20, strand = "+")
b9 <- data.frame(begin = 12, end = 25, strand = "+")
b10 <- data.frame(begin = 20, end = 25, strand = "+")
b11 <- data.frame(begin = 22, end = 25, strand = "+")

testResultsB <- c(isTRUE(all.equal(overlapping(a, b1), 0)),
                  isTRUE(all.equal(overlapping(a, b2), 1)),
                  isTRUE(all.equal(overlapping(a, b3), 1)),
                  isTRUE(all.equal(overlapping(a, b4), 1)),
                  isTRUE(all.equal(overlapping(a, b5), 1)),
                  isTRUE(all.equal(overlapping(a, b6), 1)),
                  isTRUE(all.equal(overlapping(a, b7), 1)),
                  isTRUE(all.equal(overlapping(a, b8), 1)),
                  isTRUE(all.equal(overlapping(a, b9), 1)),
                  isTRUE(all.equal(overlapping(a, b10), 1)),
                  isTRUE(all.equal(overlapping(a, b11), 0))
)

if(!isTRUE(all(testResultsB))){
  stop(paste("failed tests for function overlapping, test No.",which(testResultsB == FALSE)))
}

# test function exon2intronsBed

exons <- data.frame(scaffold = rep("sc1", 3),
                    start = c(2, 7, 12), 
                    end = c(5, 10, 20), 
                    name = c("exon1", "exon2", "exon3"), 
                    gene = rep("gene1", 3), stringsAsFactors = FALSE)
introns <- data.frame(scaffold = rep("sc1", 2),
                      intronStart = c(5, 10),
                      intronEnd = c(7, 12),
                      gene = rep("gene1", 2), stringsAsFactors = FALSE)
if(!isTRUE(all.equal(introns, exon2intronsBed(exons)))){
  stop(paste("failed tests for function exon2intronsBed"))
}

# test function gene to intergenic

genes <- data.frame(scaffold = rep("sc1", 3),
                    start = c(2, 7, 12),
                    end = c(5, 10, 20),
                    name = c("gene1", "gene2", "gene3"), stringsAsFactors = FALSE)
scaffoldLengths <- data.frame(scaffold = "sc1", length = 30)

intergenic <- data.frame(scaffold = rep("sc1", 4),
                         interStart = c(0, 5, 10, 20),
                         interEnd = c(2, 7, 12, 30), stringsAsFactors = FALSE)

if(!isTRUE(all.equal(intergenic, gene2intergenicBed(genes, scaffoldLengths)))){
  stop(paste("failed tests for function exon2intergenicBed"))
}

