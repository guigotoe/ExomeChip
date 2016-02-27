dataFileI<-commandArgs()[6]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomechipz.imiss'#
dataFileL<-commandArgs()[7]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomechipz.lmiss'#
resultFile<-commandArgs()[8]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/missingz.jpeg'#

#Missing
di=read.table(dataFileI, header=T,as.is=T)
dl=read.table(dataFileL, header=T,as.is=T)

if (length(which(di$MISS_PHENO=='Y'))<=1) {
  cat(paste("There are not Missing phenotype in MISS_PHENO of file ",dataFileI,"; Will not generate the figure\n",sp=""))
#} else {
  indy=subset(di, di$MISS_PHENO=='Y')
  jpeg('y'+resultFile, res=300, width=10, height=5, units="in")
  #par(mfcol=c(1,2), pty="m")
  hist(indy[,6], main="SNPs missing by Individuals", ylab="Freq", xlab="Missing proportion", breaks =100)
  result<-dev.off()
  cat(paste('y'+resultFile," was successfully generated.\n",sep=""))

ind=subset(di, di$MISS_PHENO=='N')
snp=subset(dl, dl$N_MISS>0)
jpeg(resultFile, res=300, width=10, height=5, units="in")
par(mfcol=c(1,2), pty="m")
hist(ind[,6], main="SNPs missing by Individuals", ylab="Freq", xlab="Missing proportion", breaks =100)
hist(snp[,3], main="SNPs missed", ylab="Freq", xlab="Number of SNPs missed", breaks=100)
result<-dev.off()
cat(paste(resultFile," was successfully generated.\n",sep=""))