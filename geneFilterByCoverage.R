# filter genes with low coverage
# and make a report with the distributions of coverage over genes

plotCovHist <- function(abr, ...){
  # load average per nucleotide read coverage over genes
  covGF <- paste("~/data/IES/analysis/coverage/", abr,".gene.cov", sep = "")
  b <- read.table(covGF, stringsAsFactors = FALSE)
  M <- max(b$V2)
  qs <- quantile(b$V2, probs = c(0.1, 0.9), names = FALSE)
  hb <- hist(log10(b$V2), breaks = 1000, xlab = "log coverage", ylab = "counts", ...)
  abline(v=log10(c(15, qs)), col = c("red", "blue", "blue"))
  filtGenesI <- which(b$V2>=15 | (b$V2 <qs[2] & b$V2 > qs[1]))
  b$V1[filtGenesI]
}

#plot

pdf(file="~/data/IES/analysis/figures/coverageOverGenes.pdf")
par(las = 1, bty = "l")
par(mfrow = c(2,2))
gppr <- plotCovHist("ppr", main = expression(italic("P. primaurelia")),   xlim = c(0, 3))
gpbi <- plotCovHist("pbi", main = expression(italic("P. biaurelia")),     xlim = c(0, 3))
gpte <- plotCovHist("pte", main = expression(italic("P. tetraurelia")),   xlim = c(0, 3))
gppe <- plotCovHist("ppe", main = expression(italic("P. pentaurelia")),   xlim = c(0, 3))
gpse <- plotCovHist("pse", main = expression(italic("P. sexaurelia")),    xlim = c(0, 3))
gpoc <- plotCovHist("poc", main = expression(italic("P. octaurelia")),    xlim = c(0, 3))
gptr <- plotCovHist("ptr", main = expression(italic("P. tredecaurelia")), xlim = c(0, 3))
gpso <- plotCovHist("pso", main = expression(italic("P. sonneborni")),    xlim = c(0, 3))
gpca <- plotCovHist("pca", main = expression(italic("P. caudatum")),      xlim = c(0, 3))


write.table(gppr, file = "~/data/IES/analysis/tables/gppr.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpbi, file = "~/data/IES/analysis/tables/gpbi.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpte, file = "~/data/IES/analysis/tables/gpte.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gppe, file = "~/data/IES/analysis/tables/gppe.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpse, file = "~/data/IES/analysis/tables/gpse.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpoc, file = "~/data/IES/analysis/tables/gpoc.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gptr, file = "~/data/IES/analysis/tables/gptr.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpso, file = "~/data/IES/analysis/tables/gpso.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
write.table(gpca, file = "~/data/IES/analysis/tables/gpca.filt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
dev.off()