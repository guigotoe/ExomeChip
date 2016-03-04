script.dir <- "/home/torres/Documents/Projects/ExomeChip/bin"
sv.dir <- '/home/torres/Documents/Projects/ExomeChip/results/03_16_Assoc/singleVariant'
source(paste(script.dir,"assoc_supp.R",sep="/"))


names <- read.table(paste(script.dir,'ExomeChip_SNPsInfo.txt',sep="/"),header=T,as.is=T,na.strings = ".")[,c(1,2,3,9)]

assoc <- read.table(paste(sv.dir,"Ageing_ExomeChip_QCed.assoc",sep="/"),header=T,as.is=T)
assoc <- assoc[complete.cases(assoc[,9]),]   # Removing rows containing NAs in pvalue
assoc$Bonf <- p.adjust(assoc$"P",method="bonferroni")
assoc <- merge(assoc,names,by.x="SNP",by.y="Chr.position")

loco <- read.table(paste(sv.dir,"gcta/data_loco.loco.mlma",sep="/"),header=T,as.is=T)
loco <- loco[complete.cases(loco[,9]),]   # Removing rows containing NAs in pvalue
loco <- merge(loco,names,by.x="SNP",by.y="Chr.position")

loco_covar <- read.table(paste(sv.dir,"gcta/dataAll_covar_loco.loco.mlma",sep="/"),header=T,as.is=T)
loco_covar <- loco_covar[complete.cases(loco_covar[,9]),]   # Removing rows containing NAs in pvalue
loco_covar <- merge(loco_covar,names,by.x="SNP",by.y="Chr.position")

snps <- read.table(paste(sv.dir,"tocheck-All-loco_covar_sig_clusterplot.txt.scores",sep="/"),header=F,as.is=T)
colnames(snps) <- c("id","score")
snps2  <- snps[snps$id==1]

k <- merge(snps,loco_covar,by.x="id",by.y="exmID",all.x=T)
k <- subset(k,k$score==1)
write.table(k,paste(sv.dir,'signf_clusterp.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")

## manhattan plot ##
require(qqman)
df <- assoc[,c(13,2,3,9)]

colnames(df) <- c("SNP","CHR","BP","P")
df[df<1e-16] <- 1e-16
manhattan(df,main = "Ageing ExomeChip Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df$P, main = "Q-Q plot of ExomeChip p-values")
tocheck <- assoc[assoc["P"]<=1e-5,]
write.table(tocheck[12],paste(sv.dir,'tocheck_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
tocheck2 <- assoc[assoc["P"]<=1e-8,]
write.table(tocheck2[12],paste(sv.dir,'signf_tocheck_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")





## MLMA-loco

df_mlma_l <- loco[,c("SNP","Chr","bp","p")]
colnames(df_mlma_l) <- c("SNP","CHR","BP","P")
df_mlma_l[df_mlma_l<1e-16] <- 1e-16
manhattan(df_mlma_l,main = "Ageing ExomeChip MLMA-loco Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df_mlma_l$P, main = "Q-Q plot of Ageing ExomeChip MLMA-loco  p-values")
#
tocheck_loco <- loco[loco["p"]<=1e-5,]
write.table(tocheck_loco["exmID"],paste(sv.dir,'tocheck-loco_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
tocheck_loco_sig <- loco[loco["p"]<=1e-8,]
write.table(tocheck_loco_sig["exmID"],paste(sv.dir,'tocheck-loco_sig_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")

## MLMA-loco-covar

df_loco_covar <- loco_covar[,c("SNP","Chr","bp","p")]
colnames(df_loco_covar) <- c("SNP","CHR","BP","P")
df_loco_covar[df_loco_covar<1e-16] <- 1e-16
manhattan(df_loco_covar,main = "Ageing ExomeChip Allcovar-MLMA-loco Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df_loco_covar$P, main = "Q-Q plot of Ageing ExomeChip covar-MLMA-loco  p-values")
#
tocheck_loco_covar <- loco_covar[loco_covar["p"]<=1e-5,]
write.table(tocheck_loco_covar["exmID"],paste(sv.dir,'tocheck-All-loco_covar_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
tocheck_loco_covar_sig <- loco_covar[loco_covar["p"]<=1e-3,]
write.table(tocheck_loco_covar_sig["exmID"],paste(sv.dir,'tocheck-All-loco_covar_sig_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
write.table(tocheck_loco_covar_sig,paste(sv.dir,'All-loco_covar_sig_SNPs.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")

ls

###
tocheck <- assoc[assoc["P"]<=0.01,]
write.table(tocheck[12],paste(sv.dir,'tocheck_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")

#
scores <- read.table('/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/ClusPlot_filtering/tocheck_clusterplot.txt.scores',header=F,as.is=T)
no <- scores[scores[2]!=1,1]
nochrpos <- names[names$'exmID'%in%no,3]
write.table(nochrpos,paste(sv.dir,'fail_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
###
write.table(assoc,paste(sv.dir,'assoc.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(assoc_signif["exmID"],paste(sv.dir,'exmId_sig.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")


## immuno chip
im.dir <- '/home/torres/Documents/Projects/ImmunoChip/results/03_16_assoc'
loco_im <- read.table(paste(im.dir,"data_covar_loco.loco.mlma",sep="/"),header=T,as.is=T)
loco_im <- loco_im[complete.cases(loco_im[,"p"]),]   # Removing rows containing NAs in pvalue
loco_im <- merge(loco_im,names,by.x="SNP",by.y="RsID",all.x=T)

df <- loco_im[,c("SNP","Chr","bp","p")]
colnames(df) <- c("SNP","CHR","BP","P")
df[df<1e-16] <- 1e-16
manhattan(df,main = "Ageing ImmunoChip Manhattan Plot", cex=0.5,cex.axis=0.8,ymax=18)
qq(df$P, main = "Q-Q plot of ExomeChip p-values")
tocheck <- loco_im[loco_im["p"]<=1e-3,]
write.table(tocheck["SNP"],paste(im.dir,'tocheck-immuno_sig_clusterplot.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
write.table(tocheck,paste(im.dir,'immuno-loco_covar_sig_SNPs.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")

ja.dir <- "/home/torres/Documents/Projects/ExomeChip/docs/janina/"
hla.dqb1 <- read.table(paste(ja.dir,"HLA-DQB1.csv",sep="/"),header=T,as.is=T,sep="\t")
loco_im_hla.dq <- subset(loco_im,loco_im$SNP %in% hla.dqb1$Variant.ID)
write.table(loco_im_hla.dq,paste(im.dir,'hla-dqa1_immuno-loco.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")

hla.dr <- read.table(paste(ja.dir,"HLA-DR.csv",sep="/"),header=T,as.is=T,sep="\t",na.strings = "",quote = "\"")
loco_im_hla.dr <- subset(loco_im,loco_im$SNP %in% hla.dr$Variant.ID)
write.table(loco_im_hla.dr,paste(im.dir,'hla-dr_immuno-loco.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")

hla.ensmbl <- read.table(paste(ja.dir,"HLA-SNPs Ensembl",sep="/"),header=T,as.is=T,sep="\t",na.strings = "",quote = "\"")
hla.ensmbl <- hla.ensmbl[complete.cases(hla.ensmbl[,1]),] 
loco_im_hla.e <- subset(loco_im,loco_im$SNP %in% hla.ensmbl$SNP)
write.table(loco_im_hla.e,paste(im.dir,'hla-ensmbl_immuno-loco.txt',sep="/"),row.names=F,col.names=T,quote=F,sep="\t")


