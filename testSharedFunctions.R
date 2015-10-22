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
