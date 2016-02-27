path <- '/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/data_qced_Test'
chid <- '/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/rawdata/New_snp_ids.txt'
f1<-'/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/rawdata/data.bim'#commandArgs()[7]#
f2<-'/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/data_qced_Test/Ageing_ExomeChip_QCed.fam'#commandArgs()[7]#
st <- '/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/stratification/data.hapmap_snps.flip.pca.evec'
snp <- '/home/torres/Documents/Projects/ExomeChip/results/02_16_QC/data_qced_Test/ExomeChip_LociNames.csv'
snp_extra <- '/home/torres/Documents/Projects/ExomeChip/docs/HumanCoreExome-12v1-1_A.annotation.txt'

## globals ##
b <- read.table(f2,header=F,as.is=T) #FamilyID-IndividualID-PaternalID-MaternalID-Sex(1=M,2=F,other)-Phenotype
c <- read.table(st,header=F,as.is=T) 
exm <- read.table(f1,header=F,as.is=T)
colnames(exm) <- c("chr","ID","unk","Pos","Maj","Min")
s <- read.table(snp,header=T,as.is=T,na.strings = ".")
se <- read.table(snp_extra,header=T,as.is=T,sep = "\t",na.strings = "")
chpid <-  read.table(chid,header=F,as.is=T)

snp_info <- merge(s,se,by.x="exmID",by.y="Name",all.x=T)
snp_info <- snp_info[,c(1,4,3,5,6,7,8,9,10,11,12,2)]
write.table(snp_info,paste(path,"ExomeChip_SNPsInfo.txt",sep="/"),quote=F,sep="\t",row.names=F)

#snp_info <- merge(s,exm,by="Chr.position")
##

## Defining Covariants gender and stratification PCA eigenvec.##
bx <- b[,c(1,2,5)]
colnames(bx) <- c("FID","IID","COV1")
cx <- c[,c(1,2,3,4)]
colnames(cx) <- c("IID","COV21","COV3","COV4")
cov <- merge(bx,cx, by = "IID")
cov <- cov[,c(2,1,3,4,5,6)]
write.table(cov,'/home/torres/Documents/Projects/ExomeChip/results/28_01_16/data_qced_Test/covariates.txt',
            row.names=F,col.names=F,quote=F,sep="\t")
##

## Get relational table Genes - SNPs; retrieve infor from ensembl ##
require(biomaRt)
listMarts(host="www.ensembl.org")
ensembl = useMart(biomart="ENSEMBL_MART_SNP",host="www.ensembl.org")
#listDatasets(ensembl)
snp_mart = useDataset("hsapiens_snp",mart=ensembl)
#attributes = listAttributes(snp_mart)
#filters = listFilters(snp_mart)
genes_snp1 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[1:61192,2],mart=snp_mart)
genes_snp2.1 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[61193:91789,2],mart=snp_mart)
genes_snp2.2 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[91790:122385,2],mart=snp_mart)
genes_snp3 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[122386:183578,2],mart=snp_mart)
genes_snp4.1 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[183579:214175,2],mart=snp_mart)
genes_snp4.2 <- getBM(attributes=c("refsnp_id","ensembl_gene_stable_id"),filters=c("snp_filter"),values=snp_info[214176:244770,2],mart=snp_mart)


genes_snp <- rbind(genes_snp1,genes_snp2.1,genes_snp2.2,genes_snp3,genes_snp4.1,genes_snp4.2)
write.table(genes_snp,paste(path,"All_SNPs_geneVariants.txt",sep="/"),quote=F,sep="\t",row.names=F)




ind <- which(with(genes_snp,ensembl_gene_stable_id == ""))
genes_snp_filtred <- genes_snp[ -ind, ]
write.table(genes_snp_filtred,paste(path,"SNPs_genes.txt",sep="/"),quote=F,sep="\t",row.names=F)
write.table(snp_info[,c(1,4)],paste(path,"ChrPos_RS_IDs.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
write.table(unique(genes_snp_filtred[1]),paste(path,"snps2assoc.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
write.table(snp_info[duplicated(snp_info[1]),1],paste(path,"snps_duplicated.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
head(genes_snp)



snp_info_all <- merge(snp_info,genes_snp,by.x="RsID",by.y="refsnp_id",all.x=T)
snp_info_all <- snp_info_all[,c(2,1,6,8,7,10,9,13,11,12,4,5,6,3)]

                                1,2,3,4,5,6,7,8,9,10,11,12,13
##

## Using SKAT for gene based asociation test ##
require(SKAT)
bed <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/data_qced_Test/plink.bed'#commandArgs()[7]#
bim <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/data_qced_Test/plink.bim'#commandArgs()[8]#
fam <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/data_qced_Test/plink.fam'#commandArgs()[9]#
cov <- '/home/torres/Documents/Projects/ExomeChip/results/28_01_16/data_qced_Test/plink.cov'#commandArgs()[9]#
setID <- 

ssd <- Generate_SSD_SetID(bed,bim,fam)
