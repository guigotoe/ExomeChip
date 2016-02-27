#fileName<-commandArgs()[6]
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_ALL.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_Case.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_Control.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCDEBONN.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCDEMICK.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCGCON.hwe"
fileName <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCPOPGEN.hwe"

#resultFile<-paste0(fileName,".HWE.jpg")
fileContent<-read.table(fileName, header=T,as.is=T)
head(fileContent)

p.hwe.cont <- fileContent$"P"[fileContent$"TEST"=="UNAFF"]
p.hwe.case <- fileContent$"P"[fileContent$"TEST"=="AFF"]

par(mfrow=c(1,2))
x <- sort(-log10(ppoints(p.hwe.cont))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.cont))
plot(x,y,ylab="Observed",xlab="Expected",main="Controls")
abline(a=0,b=1,lty=2)

x <- sort(-log10(ppoints(p.hwe.case))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.case))
plot(x,y,ylab="Observed",xlab="Expected",main="Cases")
abline(a=0,b=1,lty=2)

par(mfrow=c(1,1))




data.sub1 <- subset(fileContent,fileContent$"O.HET."<0.05)
data.sub2 <- subset(fileContent,fileContent$"P"<0.00005)


p.hwe.cont <- fileContent$"P"[fileContent$"TEST"=="UNAFF"]
#pbf.hwe.cont <- p.adjust(p.hwe.cont,method="bonferroni",n=length(p.hwe.cont))
p.hwe.case <- fileContent$"P"[fileContent$"TEST"=="AFF"]
#pbf.hwe.case <- p.adjust(p.hwe.case,method="bonferroni",n=length(p.hwe.case))
par(mfrow=c(1,2))
# x <- sort(-log10p(points(p.hwe.cont))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.cont))
plot(x,y,ylab="Observed",xlab="Expected",main="Controls")
abline(a=0,b=1,lty=2)

x <- sort(-log10(ppoints(p.hwe.case))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.case))
plot(x,y,ylab="Observed",xlab="Expected",main="Cases")
abline(a=0,b=1,lty=2)

#library(geneplotter)
colors  <- densCols(data.sub1$"O.HET.")
y <- sort(data.sub1[,7]) #p.hwe pvalue list
par(mfrow=0)

jpeg(resultFile, res=300, width=5, height=5, units="in")
hist(fileContent[,8],main="HWE Distribution",xlab="",las=1)
hist(fileContent[,9],main="HWE P-value",xlab="",las=1)
hist(fileContent[,7],main="Heterosity rate",xlab="",las=1)
hist(data.sub1[,7],main="Heterosity rate < 0.05",xlab="",las=1)
plot(x)


dev.off()
cat(paste(resultFile," was successfully generated.\n",sep=""))
