library(car)
library(ggplot2)
dataFile<-commandArgs()[6]#'/home/torres/Documents/Projects/ExomeChip/data/zcalls/pca_noexclusions/exomez_short_m0.pca.evec'#commandArgs()[6]

data <- read.table(dataFile, header=F)
colnames(data)<-c("id","I","II","III","IV","V","VI","VII","VIII","IX","X","Type")
levels(data$Type)[levels(data$Type)==3] <- "CEU"
levels(data$Type)[levels(data$Type)==4] <- "CHB"
levels(data$Type)[levels(data$Type)==5] <- "JPT"
levels(data$Type)[levels(data$Type)==6] <- "YRI"
data <- data.frame(data)
data.sub <- subset(data,!data$Type == "YRI")
data.sub2 <- subset(data,data$Type %in% c("CEU","Case","Control"))
data.sub3 <- subset(data, data$Type %in% c("Control","Case"))

png("eigenstratPCA1.png", res=300, width=10, height=5, units="in")

ggplot(data,aes(x=I,y=II,color=Type))+geom_point(shape=20)+
  xlab("PC1")+ylab("PC2")+ggtitle("Population stratification")+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
  scale_color_discrete(name="Populations")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10, face="bold"))
png("eigenstratPCA2.png", res=300, width=10, height=5, units="in")

ggplot(data.sub,aes(x=I,y=II,color=Type))+geom_point(shape=20)+
  xlab("PC1")+ylab("PC2")+ggtitle("Population stratification")+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
  scale_color_discrete(name="Populations")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10, face="bold"))

png("eigenstratPCA3.png", res=300, width=10, height=5, units="in")

ggplot(data.sub2,aes(x=I,y=II,color=Type))+geom_point(shape=20)+
  xlab("PC1")+ylab("PC2")+ggtitle("Population stratification")+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
  scale_color_discrete(name="Populations")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10, face="bold"))

png("eigenstratPCA4.png", res=300, width=10, height=5, units="in")

ggplot(data.sub3,aes(x=I,y=II,color=Type))+geom_point(shape=20)+
  xlab("PC1")+ylab("PC2")+ggtitle("Population stratification")+
  geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
  scale_color_discrete(name="Populations")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10, face="bold"))

#data.sub4 <- subset(data, data$"Type" %in% c("Control","Case") & data$"I"<0)$'id'

#lapply(data.sub4, write, "failstrat.txt", append=TRUE)

