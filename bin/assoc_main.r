####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 1 March 2015
# Created: 23 February 2015
#
# This is written as part of ExomeChip analysis, 
# but could be splitted to serve different purposes.
####################################################
#script.dir <- getwd()
script.dir <- "/home/torres/Documents/Projects/ExomeChip/bin"
out.dir <- '/home/torres/Documents/Projects/ExomeChip/results/03_16_Assoc/geneBased'
path <- '/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/data_qced_test'

source(paste(script.dir,"assoc_supp.R",sep="/"))

cov <- cov_file(paste(out.dir,'data.fam',sep="/"),'/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/stratification/data.hapmap_snps.flip.pca.evec')
write.table(cov,paste(out.dir,'covariates.txt',sep="/"),row.names=F,col.names=F,quote=F,sep="\t")

bed <- paste(out.dir,'data.bed',sep="/")       #commandArgs()[7]#
bim <- paste(out.dir,'data.bim',sep="/")       #commandArgs()[8]#
fam <- paste(out.dir,'data.fam',sep="/")       #commandArgs()[9]#
cov <- paste(out.dir,'covariates.txt',sep="/") #commandArgs()[9]#

snps <- snpinfo(bim,snpid = paste(path,"ExomeChip_LociNames.csv",sep="/"),snpgene = paste(path,"SNPs_genes.txt",sep="/"))
