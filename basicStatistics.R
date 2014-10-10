
#while waiting some basic statistics
pb<-read.table("pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.ies",as.is=T)
ps1<-read.table("psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.ies",as.is=T,nrows = 24734)
ps2<-read.table("psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.ies",as.is=T,skip = 24736)

pt<-read.table("ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.ies",as.is=T)

par(bty="l",mfrow=c(2,2))
a<-hist(nchar(pb$V5),breaks=10000,xlim=c(20,200),main = expression(italic("P. biaurelia")),xlab="IES length (bp)")
identify()
hist(c(nchar(ps1$V5),nchar(ps1$V5)),breaks=10000,xlim=c(20,200),main = expression(italic("P. sexaurelia")),xlab="IES length (bp)")
hist(nchar(pt$V5),breaks=10000,xlim=c(20,200),main = expression(italic("P. tetraurelia")),xlab="IES length (bp)")

#also read files and make genebank format
