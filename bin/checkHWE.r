f1 <- commandArgs()[6]
#f1 <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/HWE/data.iqc.lmiss.hardy.hwe'#commandArgs()[6]#
all <-read.table(f1, header=T,as.is=T)
p.hwe <- subset(all, all$"TEST"=="ALL")
p.hwe["y"] <- -log10(p.hwe$"P")
p.hwe <- p.hwe[order(p.hwe$y),]
p.hwe["x"] <- sort(-log10(ppoints(p.hwe$"P")))#expected pvalue list

p.hwe.cont <- subset(all,all$"TEST"=="UNAFF")
p.hwe.cont["BONF"] <- p.adjust(p.hwe.cont$"P",method="bonferroni",n=NROW(p.hwe.cont))
p.hwe.cont["y"] <- -log10(p.hwe.cont$"P")
p.hwe.cont <- p.hwe.cont[order(p.hwe.cont$y),]
p.hwe.cont["x"] <- sort(-log10(ppoints(p.hwe.cont$"P")))#expected pvalue list

p.hwe.case <- subset(all,all$"TEST"=="AFF")
p.hwe.case["BONF"] <- p.adjust(p.hwe.case$"P",method="bonferroni",n=NROW(p.hwe.case))
p.hwe.case["y"] <- -log10(p.hwe.case$"P")
p.hwe.case <- p.hwe.case[order(p.hwe.case$y),]
p.hwe.case["x"] <- sort(-log10(ppoints(p.hwe.case$"P"))) #expected pvalue list
  
jpeg(sprintf("%s_qqplot.jpg",f1), res=300, width=10, height=5, units="in")
par(mfrow=c(1,2))
plot(p.hwe.cont$"O.HET.",p.hwe.cont$"E.HET.",ylab="Observed",xlab="Expected",
     main="HWE-test QQplot - Controls",pch=20,col=ifelse(p.hwe.cont$BONF<=0.05,"red","black"))
abline(a=0,b=1,lty=2,col='red')
plot(p.hwe.case$"O.HET.",p.hwe.case$"E.HET.",ylab="Observed",xlab="Expected",
     main="HWE-test QQplot - Cases",pch=20,col=ifelse(p.hwe.case$BONF<=0.05,"red","black"))
abline(a=0,b=1,lty=2,col='red')
dev.off()
jpeg(sprintf("%s_bonf_hist.jpg",f1), res=300, width=10, height=5, units="in")
par(mfrow=c(1,1))
bonf <- c(p.hwe.cont$BONF,p.hwe.case$BONF)
hist(bonf,main="HWE test corrected P value (Bonf.)",xlab="P value")
dev.off()
fails  <- unique(c(subset(p.hwe.cont,p.hwe.cont$BONF<=0.05)$SNP,subset(p.hwe.case,p.hwe.case$BONF<=0.05)$SNP))
write.table(fails,file=sprintf("%s_fails.txt",f1),quote=FALSE,sep="\n",na="NA",row.names=FALSE)

cat(paste("Graphs and fails file were successfully generated.\n",sep=""))  
