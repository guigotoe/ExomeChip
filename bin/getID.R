####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 7 March 2015
# Created: 23 February 2015
#
# Supplementary functions - getID.R.
# For help use: getID.R -h
# This is written as part of ExomeChip analysis, 
# but could be splitted to serve different purposes.
####################################################
# 1 = "exmID", 2 = "RsID", 3 = "Chr.position", 4 = "RsID_Chr.pos", 5 = "Chr",
# 6 = "MapInfo", 7 ="Alleles.y",8 = "Transcript.s.", 8 = "Gene.s.", 9= "In.exon",
# 10 = "Mutation.s.", 11 = "Alleles.x")

getID <- function(query,col,request,ref,header=FALSE){
  colnames(query) <- col
  info <- merge(query,ref,by=col)
  write.table(info[request],paste(qname,'newid.txt',sep="."),col.names=header,row.names=F,sep="\t",na="NA")
}

script.dir <- "/home/torres/Documents/Projects/ExomeChip/bin"
names <- read.table(paste(script.dir,'ExomeChip_SNPsInfo.txt',sep="/"),header=T,as.is=T) # IDs reference file
qname <- commandArgs()[6]#'./Var_chr-position.txt'  # query file name
inp <- commandArgs()[7]#"2,1,3" # Info. requested (list of numbers separated by comma)
query <- read.table(qname,header=F,as.is=T)  # query list of IDs

if (length(grep("^[(exm)]",query))==1){col="exmID"}else if(length(grep("chr[[:digit:]]+:",query))==1){
  col="Chr.position"}else if(length(grep("^[(rs)]",query))==1){col="RsID"}
request=as.numeric(unlist(strsplit(inp,',')))
getID(query,col,request,ref=names)


