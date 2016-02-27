f1 <- commandArgs()[6]
#f1<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez.sex.imisshet.nameupd.strat.ibd.lmiss'#
#f2<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez.sex.imisshet.nameupd.strat.ibd.imiss'#

l=read.table(f1, header=T,as.is=T)
#i=read.table(f2, header=T,as.is=T)
x <- sort(l$"F_MISS")#log10(l$"F_MISS"))
x <- x[x>=0.02]

jpeg("SNPmissing.jpg", res=300, width=10, height=5, units="in")
hist(x, breaks =100,main="SNPs missing rate > 0.02",ylab="Number of SNPs",xlab="Fraction of missing data")
abline(v=0.03,col="RED",lty=2)

cat(paste("The graph was successfully generated.\n",sep=""))
