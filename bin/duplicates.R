b <- commandArgs()[6]
bim <- read.table(b,header=F,as.is=T)
dup <- bim[duplicated(bim[2]),2]
write.table(dup,"duplicates.txt",quote=F,sep="\t",row.names=F,col.names=F)