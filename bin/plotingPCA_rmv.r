library(car)
dataFile<-commandArgs()[6]#
dataFile2<-commandArgs()[7]#
rmfile <- commandArgs()[8]#
#dataFile<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez_short.pca.evec'#commandArgs()[6]
#rmfile <- '/home/torres/Documents/Projects/ExomeChip/data/zcalls/pca_noexclusions/failstrat.txt'#commandArgs()[7]
#dataFile2<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/pca_noexclusions/exomez_short_m0.pca.evec'#commandArgs()[6]

data2 <- read.table(dataFile2,header=F)
data <- read.table(dataFile, header=F)
rm <- read.table(rmfile, header=F)

colnames(data)<-c("id","I","II","III","IV","V","VI","VII","VIII","IX","X","Type")
levels(data$Type)[levels(data$Type)==3] <- "CEU"
levels(data$Type)[levels(data$Type)==4] <- "CHB"
levels(data$Type)[levels(data$Type)==5] <- "JPT"
levels(data$Type)[levels(data$Type)==6] <- "YRI"
data <- data.frame(data)


data$"batch"<-NA
data$"batch"[grep("AGI",data$"id")] <- "AGE"
data$"batch"[grep("HCDEBONN",data$"id")]<-"HCDEBONN"
data$"batch"[grep("HCDEKORA2",data$"id")]<-"HCDEKORA2"
data$"batch"[grep("HCDEMICK",data$"id")]<-"HCDEMICK"
data$"batch"[grep("HCGCON",data$"id")]<-"HCGCON"
data$"batch"[grep("HCPOPGEN",data$"id")]<-"HCPOPGEN"


jpeg("eigenstratPCA1.1.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data,
            legend.title="Samples",legend.coord="topright")
result<-dev.off()

jpeg("eigenstratPCA1_batches.1.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | batch, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data,
            legend.title="Batches",legend.coord="topright")
result<-dev.off()

data.sub <- subset(data,!(data$"id" %in% rm$"V1"))


jpeg("eigenstratPCA1.2.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data.sub,
            legend.title="Samples",legend.coord="topright")
result<-dev.off()

rmlst2 <- subset(data, data$"Type" %in% c("Control","Case") & data$"I">0.03)$"id"
data.sub2 <- subset(data,!data$"id" %in% rmlst2)


jpeg("eigenstratPCA1.3.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data.sub2,
            legend.title="Samples",legend.coord="topright")
result<-dev.off()

####

colnames(data2)<-c("id","I","II","III","IV","V","VI","VII","VIII","IX","X","Type")
levels(data2$Type)[levels(data2$Type)==3] <- "CEU"
levels(data2$Type)[levels(data2$Type)==4] <- "CHB"
levels(data2$Type)[levels(data2$Type)==5] <- "JPT"
levels(data2$Type)[levels(data2$Type)==6] <- "YRI"
data2 <- data.frame(data2)
data2.sub1 <- subset(data2,data2$"Type" %in% c("CEU","Case","Control") & data2$"I">-0.01)
data2.sub3 <- subset(data2, data2$"Type" %in% c("CEU","Control","Case") & data2$"I">-0.04 & !(data2$"id" %in% rm$"V1"))
data2.sub4 <- subset(data2.sub3, !(data2.sub3$"id" %in% rmlst2))

jpeg("eigenstratPCA2.1.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data2,
            legend.title="Samples",legend.coord="topleft")
result<-dev.off()

jpeg("eigenstratPCA2.2.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data2.sub1,
            legend.title="Samples",legend.coord="topleft")
result<-dev.off()

jpeg("eigenstratPCA2.3.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data2.sub3,
            legend.title="Samples",legend.coord="topleft")
result<-dev.off()

jpeg("eigenstratPCA2.4.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data2.sub4,
            legend.title="Samples",legend.coord="topleft")
result<-dev.off()
rm2<- lapply(rmlst2, as.character)
rm1<- lapply(rm$"V1", as.character)
rmlist <- unique(c(rm1,rm2))

write.table(rmlist, "failstrat.txt", sep="\n",quote=FALSE,row.names=FALSE,col.names=FALSE) 
