fmap<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomechipz.map'#commandArgs()[6]#
fped<-'/home/torres/Documents/Projects/ExomeChip/data/zcalls/exomechipz.ped'#commandArgs()[6]#

map <- read.delim(fmap,header=F,as.is=T)
head(map)
nrow(map)

#qqplot
# we're using the P points command to estimate the quantiles of a uniform distribution, which
# is the distribution that the P values should have under the null. 
#x<-sort(-log10(ppoints(p.hwe)))#p.hwe pvalue list
#y<-sort(-log10(p.hwe)) #p.hwe pvalue list
#plot(x,y,ylab="Observed",xlab="expected",main='All')
#abline(a=0,b=1,lty=2)