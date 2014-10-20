
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


#also extract sequencies and make blastdb and all against all alignment and run them on silix to see how similar they are
#length distribution
#gc% vs length
#?assemble

#make fasta files with proteins in each group from silix output
run t_coffee to align them
t_coffee foo.seq -mode accurate

perl basicStatistics.pl ~/data/IES_data/pbiaurelia/PBI.ies >  ~/data/IES_data/pbiaurelia/Pbi.length.gc.tab &
perl basicStatistics.pl ~/data/IES_data/psexaurelia/PSEX.ies >  ~/data/IES_data/psexaurelia/Pse.length.gc.tab &
perl basicStatistics.pl ~/data/IES_data/ptetraurelia/PTET.ies >  ~/data/IES_data/ptetraurelia/Pte.length.gc.tab

#length and GC% are positively correlated.

Pbi<-read.table("pbiaurelia/Pbi.length.gc.tab")
Pte<-read.table("ptetraurelia/Pte.length.gc.tab")
Pse<-read.table("psexaurelia/Pse.length.gc.tab")

par(bty="l",mfrow=c(2,2))
plot(Pbi$V2,Pbi$V3,pch=19,cex=0.5,main=expression(italic("P. biaurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pbi$V3 ~ Pbi$V2)
abline(fit,col="red")

plot(Pte$V2,Pte$V3,pch=19,cex=0.5,main=expression(italic("P. tetraurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pte$V3 ~ Pte$V2)
abline(fit,col="red")

plot(Pse$V2,Pse$V3,pch=19,cex=0.5,main=expression(italic("P. sexaurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pse$V3 ~ Pse$V2)
abline(fit,col="red")

# I exclude the small sizes because of noise and see the same pattern
index<-which(Pbi$V2>1000)
Pbi<-Pbi[index,]
index<-which(Pte$V2>1000)
Pte<-Pte[index,]
index<-which(Pse$V2>1000)
Pse<-Pse[index,]

par(bty="l",mfrow=c(2,2))
plot(Pbi$V2,Pbi$V3,pch=19,cex=0.5,main=expression(italic("P. biaurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pbi$V3 ~ Pbi$V2)
abline(fit,col="red")

plot(Pte$V2,Pte$V3,pch=19,cex=0.5,main=expression(italic("P. tetraurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pte$V3 ~ Pte$V2)
abline(fit,col="red")

plot(Pse$V2,Pse$V3,pch=19,cex=0.5,main=expression(italic("P. sexaurelia")),xlab="IES length",ylab="GC%")
fit<-lm(Pse$V3 ~ Pse$V2)
abline(fit,col="red")

#If small IES are older they seem to start as high in GC and get closser to the genome average:
./describe_seq.exe ~/data/IES_data/pbiaurelia/biaurelia_V1-4_assembly_v1.fasta
Using file:/Users/diamantis/data/IES_data/pbiaurelia/biaurelia_V1-4_assembly_v1.fasta

scaffold_0001
Number of residues
A	27913722	 36.254% A
T	27848707	 36.169% T
C	 9666917	 12.555% C
G	 9702736	 12.602% G
N	 1849232	  2.402% N
other	   14166	  0.018% other
G+C	19369653	 25.157% G+C
A+T	55762429	 72.423% A+T
CpG	  571007	  0.742% CpG
TpA	10444639	 13.565% TpA
Sequence length = 76995480

umr55558-etu-duret:c diamantis$ ./describe_seq.exe ~/data/IES_data/ptetraurelia/ptetraurelia_mac_51.fa
Using file:/Users/diamantis/data/IES_data/ptetraurelia/ptetraurelia_mac_51.fa

scaffold51_9
Number of residues
A	25809866	 35.793% A
T	25780491	 35.752% T
C	10038176	 13.921% C
G	10063901	 13.957% G
N	  411899	  0.571% N
other	    4176	  0.006% other
G+C	20102077	 27.878% G+C
A+T	51590357	 71.545% A+T
CpG	  592642	  0.822% CpG
TpA	 9320036	 12.925% TpA
Sequence length = 72108509

umr55558-etu-duret:c diamantis$ ./describe_seq.exe ~/data/IES_data/psexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta
Using file:/Users/diamantis/data/IES_data/psexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta

scaffold_001
Number of residues
A	25503116	 37.491% A
T	25497567	 37.483% T
C	 8088230	 11.890% C
G	 8080939	 11.879% G
N	  851962	  1.252% N
other	    3276	  0.005% other
G+C	16169169	 23.769% G+C
A+T	51000683	 74.973% A+T
CpG	  393833	  0.579% CpG
TpA	 9526347	 14.004% TpA
Sequence length = 68025090