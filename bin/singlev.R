script.dir <- "/home/torres/Documents/Projects/ExomeChip/bin"
sv.dir <- '/home/torres/Documents/Projects/ExomeChip/results/03_16_Assoc/singleVariant'

names <- read.table(paste(script.dir,'ExomeChip_SNPsInfo.txt',sep="/"),header=T,as.is=T,na.strings = ".")[,c(1,2,3,9)]
assoc <- read.table(paste(sv.dir,"Ageing_ExomeChip_QCed.assoc",sep="/"),header=T,as.is=T)
assoc <- assoc[complete.cases(assoc[,9]),]   # Removing rows containing NAs in pvalue
assoc$Bonf <- p.adjust(assoc$"P",method="bonferroni")
assoc <- merge(assoc,names,by.x="SNP",by.y="Chr.position")

## manhattan plot ##
require(qqman)
df <- assoc[,c(13,2,3,9)]
colnames(df) <- c("SNP","CHR","BP","P")
df[df<1e-16] <- 1e-16
manhattan(df,main = "Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df$P, main = "Q-Q plot of ExomeChip p-values")


assoc2 <- read.table('/home/torres/Documents/Projects/AffyChip/assoc_alle_individuen_60_neu.assoc',header=T,as.is=T)
assoc2 <- assoc2[complete.cases(assoc2[,9]),]   # Removing rows containing NAs in pvalue
df2 <- assoc2[,c(2,1,3,9)]
df2[df2<1e-16] <- 1e-16
manhattan(df2,main = "Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df2$P, main = "Q-Q plot of ExomeChip p-values")


####

tocheck <- assoc[assoc["P"]<=0.01,]
write.table(tocheck[12],paste(sv.dir,'tocheck_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")


write.table(assoc,paste(sv.dir,'assoc.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(assoc_signif["exmID"],paste(sv.dir,'exmId_sig.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")

