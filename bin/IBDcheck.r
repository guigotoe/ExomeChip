library(car)
#dataFile<-commandArgs()[6]#
dataFile<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomez_sex_imisshet_nameupd_strat_IBD.genome'#commandArgs()[6]
data <- read.table(dataFile, header=F)
data[,1:10]
colnames(data)<-c("id","I","II","III","IV","V","VI","VII","VIII","IX","X","Type")
levels(data$Type)[levels(data$Type)==3] <- "CEU"
levels(data$Type)[levels(data$Type)==4] <- "CHB"
levels(data$Type)[levels(data$Type)==5] <- "JPT"
levels(data$Type)[levels(data$Type)==6] <- "YRI"
data <- data.frame(data)
data.sub <- subset(data,data$"I">-0.04)
data.sub2 <- subset(data,data$"Type" %in% c("CEU","Case","Control") & data$"I">-0.01)

jpeg("eigenstratPCA1.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data,
            legend.title="Groups")#,legend.coord="topright")
result<-dev.off()
jpeg("eigenstratPCA2.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data.sub,
            legend.title="Groups")#,legend.coord="topright")
jpeg("eigenstratPCA3.jpg", res=300, width=10, height=5, units="in")
scatterplot(II~I | Type, reg.line=FALSE, smooth=FALSE, spread=FALSE, 
            boxplots=FALSE, span=0.5, xlab="PC1", ylab="PC2", 
            main="Population stratification", cex=0.5, by.groups=TRUE, data=data.sub2,
            legend.title="Groups")#,legend.coord="topright")

data.sub3 <- subset(data, data$"Type" %in% c("Control","Case") & data$"I"<0)$'id'


write.table(data.sub3, "failstrat.txt", sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE) 
