####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 1 March 2015
# Created: 23 February 2015
#
# Supplementary functions.
# This is written as part of ExomeChip analysis, 
# but could be splitted to serve different purposes.
####################################################
script.dir <- "/home/torres/Documents/Projects/ExomeChip/bin"
## Defining Covariants gender and stratification PCA eigenvec.##
cov_file <- function(famf,pcaf){
  fam <- read.table(famf,header=F,as.is=T) #FamilyID-IndividualID-PaternalID-MaternalID-Sex(1=M,2=F,other)-Phenotype
  pca <- read.table(pcaf,header=F,as.is=T)
  bx <- fam[,c(1,2,5)]
  colnames(bx) <- c("FID","IID","COV1")
  cx <- pca[,c(1,2,3,4,5)]
  colnames(cx) <- c("IID","COV2","COV3","COV4","COV5")
  cov <- merge(bx,cx, by = "IID")
  cov <- cov[,c(2,1,3,4,5,6,7)]
  return(cov)
}
##
##
snpinfo <- function(bim,snpid,snpgene){
  exm <- read.table(bim,header=F,as.is=T)
  colnames(exm) <- c("chr","RsID","unk","Pos","Maj","Min")
  s <- read.table(snpid,header=T,as.is=T,na.strings = ".")
  merged <- merge(s,exm,by="RsID")
  sg <- read.table(snpgene,header=T,as.is=T)
  all <- merge(merged,sg,by.x="RsID",by.y="refsnp_id")
  return(all)
}
