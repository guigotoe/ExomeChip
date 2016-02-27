f1 <- '/home/torres/Documents/Projects/ExomeChip/data/calls/exomechip.frq'#commandArgs()[6]#
outf <- '/home/torres/Documents/Projects/ExomeChip/data/calls/exomechipMAF.jpeg'#commandArgs()[7]#

maf<-read.table(f1,h=T,as.is=T)
cut01<-maf[which(maf$MAF > 0.01),]
cut1<-maf[which(maf$MAF > 0.1),]
y<-sort(log10(maf$'MAF')) #p.hwe pvalue list
jpeg(outf, res=300, width=10, height=5, units="in")
plot(y,ylab='MAF',xlab=sprintf("Samples (%d)",NROW(maf)), main='MAF frequency',axes=F)
axis(2,at=c(0,-1,-2,-3,-4,-5),labels=c(1,0.1,0.01,0.001,0.0001,0.00001),tick=T)
axis(1,at=NULL)
abline(h=(log10(0.01)),lty=2,col='Red')
text(25000, -1.8,sprintf("SNP %s%% (%d)",round(100*NROW(cut01)/NROW(maf),digits=1),NROW(cut01)), col = "red")
abline(h=(log10(0.1)),lty=2,col='Blue' )
text(25000, -0.8,sprintf("SNP %s%% (%d)",round(100*NROW(cut1)/NROW(maf),digits=1),NROW(cut1)), col = "Blue")
result<-dev.off()
cat(paste(outf," was successfully generated.\n",sep=""))
