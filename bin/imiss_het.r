f1<-commandArgs()[6]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez.sex.imiss'#
f2<-commandArgs()[7]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez.sex.het'#
sdth<- commandArgs()[8]#3.2# int or float: times of sd threshold default 3.2#
pmiss<- commandArgs()[9]#0.001# float: portion of missing rates default 0.01#
outf<-commandArgs()[10]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez.sexSNPmiss.jpeg'#


imiss=read.table(f1,h=T,as.is=T)
het=read.table(f2,h=T,as.is=T)
het$hetRate = round((het$N.NM. - het$O.HOM.)/het$N.NM., 5)
x <- log10(imiss$'F_MISS')
colors  <- densCols(x,het$hetRate)
colors2  <- densCols(x,het$hetRate,colramp=colorRampPalette(c("red","white")))

i <- imiss[log10(imiss$"F_MISS")>log10(as.numeric(pmiss)),]
d <- unlist(apply(imiss,1,function(x) if(log10(as.numeric(x[6]))>log10(as.numeric(pmiss))){return(x[2])}))
e <- het[het$hetRate>=(mean(het$hetRate)+(as.numeric(sdth)*sd(het$hetRate)))| het$hetRate<=(mean(het$hetRate)-(as.numeric(sdth)*sd(het$hetRate))),]
failed <- het[het$hetRate>=(mean(het$hetRate)+(as.numeric(sdth)*sd(het$hetRate)))| het$hetRate<=(mean(het$hetRate)-(as.numeric(sdth)*sd(het$hetRate)))|het$IID%in%d,]

jpeg(outf, res=300, width=10, height=5, units="in")
plot(x,het$hetRate,col=ifelse(het$hetRate %in% e$hetRate,colors2,ifelse(x %in% log10(i$F_MISS),colors2,colors)),pch=20,xlim=c(-6,0),ylim=c(0.04,0.1),xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(1,at=c(-6,-5,-4,-3,-2,-1,0),labels=c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1))
axis(2,at=c(0.05,0.06,0.07,0.08,0.09,0.1),tick=T)
abline(h=mean(het$hetRate)-(as.numeric(sdth)*sd(het$hetRate)),col="RED",lty=2)
abline(h=mean(het$hetRate)+(as.numeric(sdth)*sd(het$hetRate)),col="RED",lty=2)
abline(v=log10(as.numeric(pmiss)), col="RED", lty=2)
#text(-1.7,0.04,sprintf("%s",round(244770*0.01)), col = "Blue")
result<-dev.off()

write.table(failed,'fail.imisshet', sep="\t", quote=F, row.names =F)
write.table(failed[,c(1,2)],"failed-imisshet-qc.txt",sep="\t", quote=F, row.names =F, col.names=F)
cat(paste(outf," and failed list was successfully generated.\n",NROW(failed)," individuals failed heterozygosity and missigness\n",sep=""))
