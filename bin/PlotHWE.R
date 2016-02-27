f1 <-commandArgs()[6]
#f1 <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/HWE/data.iqc.lmiss.indepSNP.hardy.hwe'
#f1 <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_ALL.hwe"
fileContent<-read.table(f1, header=T,as.is=T)
#head(fileContent)

cont <- subset(fileContent,fileContent$"TEST"=="UNAFF")
cont["Bonf"] <- p.adjust(cont$"P",method="bonferroni")
case <- subset(fileContent,fileContent$"TEST"=="AFF")
case["Bonf"] <- p.adjust(case$"P",method="bonferroni")
## According Yan Guo et al. 2014 - concervative p-value threshold P<1x10-5

jpeg(sprintf("%s_qqplot.jpg",f1), res=300, width=10, height=5, units="in")
par(mfrow=c(1,2))
plot(cont$"E.HET",cont$"O.HET",ylab="Observed",xlab="Expected",main="Controls",
     pch=20,col=ifelse(cont$"Bonf"<=0.05, "red", "black"))
abline(a=0,b=1,lty=2)
plot(case$"E.HET",case$"O.HET",ylab="Observed",xlab="Expected",main="Cases",
     pch=20,col=ifelse(case$"Bonf"<=0.05, "red", "black"))
abline(a=0,b=1,lty=2)
dev.off()
par(mfrow=c(1,1))

jpeg(sprintf("%s_HetRate.jpg",f1), res=300, width=10, height=5, units="in")
hist(c(case$"O.HET.",cont$"O.HET."),main="Heterosity rate",xlab="Observed heterosity",las=1)
dev.off()

jpeg(sprintf("%s_bonf_hist.jpg",f1), res=300, width=10, height=5, units="in")
hist(c(case$Bonf,cont$Bonf),main="HWE test corrected P value (Bonf.)",xlab="P value")
dev.off()

fails <- cont$SNP[cont$"Bonf"<=0.05]
write.table(fails,file="controls_HWE_fails.txt",quote=FALSE,sep="\n",row.names=FALSE,col.names = FALSE)
cat(paste("Fail list contains those HWE bonf(p-value) <0.05\nGraphs and fail-file was successfully generated.\n",sep=""))

