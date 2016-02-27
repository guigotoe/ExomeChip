#fileName<-commandArgs()[6]
Fall <- commandArgs()[6]
Fcases <- commandArgs()[7]
Fcontrol <- commandArgs()[8]
Fba <- commandArgs()[9]
Fbb <- commandArgs()[10]
Fbc <- commandArgs()[11]
Fbd <- commandArgs()[12]
Fbe <- commandArgs()[13]
#Fall <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_ALL.hwe"
#Fcases <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_Case.hwe"
#Fcontrol <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_Control.hwe"
#Fba <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCDEBONN.hwe"
#Fbb <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCDEKORA2.hwe"
#Fbc <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCDEMICK.hwe"
#Fbd <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCGCON.hwe"
#Fbe <- "/home/torres/Documents/Projects/ExomeChip/data/zcalls/eexomezIQCgeno_commonsHWE_HCPOPGEN.hwe"

#resultFile<-paste0(fileName,".HWE.jpg")
all <-read.table(Fall, header=T,as.is=T)
p.hwe.cont <- all$"P"[all$"TEST"=="UNAFF"]
p.hwe.case <- all$"P"[all$"TEST"=="AFF"]
jpeg("ALL_hwe.jpg", res=300, width=10, height=5, units="in")
par(mfrow=c(1,2))
x <- sort(-log10(ppoints(p.hwe.cont))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.cont))
plot(x,y,ylab="Observed",xlab="Expected",main="SNPs MAF>0.05 HWE-test P-value - Controls")
abline(a=0,b=1,lty=2)
x <- sort(-log10(ppoints(p.hwe.case))) #p.hwe pvalue list
y <- sort(-log10(p.hwe.case))
plot(x,y,ylab="Observed",xlab="Expected",main="SNPs MAF>0.05 HWE-test P-value Cases")
abline(a=0,b=1,lty=2)

par(mfrow=c(1,1))
p.hwe <- all$"P"[all$"TEST"=="ALL"]
jpeg("ALL_hwePvalue.jpg", res=300, width=5, height=5, units="in")
hist(p.hwe,main="HWE test P value",xlab="P value")

case <-read.table(Fcases, header=T,as.is=T)
p.hwe<- case$"P"[case$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("Cases_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="SNPs MAF>0.05 HWE-test Cases QQ-plot")
abline(a=0,b=1,lty=2)
case.out <- subset(case,case$TEST=="ALL" & case$P<=0.00001)$SNP

#all.out <- subset(all,all$TEST=="ALL" & all$P<=0.000000001)$SNP
#all40.out <- subset(all,all$TEST=="ALL" & all$P<=0.00001)$SNP
#case.out <- subset(case,case$TEST=="ALL" & case$P<=0.000000001)$SNP

control <-read.table(Fcontrol, header=T,as.is=T)
p.hwe<- control$"P"[control$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("Controls_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="SNPs MAF>0.05 HWE-test P-value Controls QQ-plot")
abline(a=0,b=1,lty=2)


ba <-read.table(Fba, header=T,as.is=T)
p.hwe<- ba$"P"[ba$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("HCDEBONN_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="HCDEBONN batch QQ-plot")
abline(a=0,b=1,lty=2)
ba.out <- subset(ba,ba$TEST=="ALL" & ba$P<=0.00001)$SNP

bb <-read.table(Fbb, header=T,as.is=T)
p.hwe<- bb$"P"[bb$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("HCDEKORA2_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="HCDEKORA2 batch QQ-plot")
abline(a=0,b=1,lty=2)
bb.out <- subset(bb,bb$TEST=="ALL" & bb$P<=0.00001)$SNP
#bb40.out <- subset(bb,bb$TEST=="ALL" & bb$P<=0.00005)$SNP

bd <-read.table(Fbd, header=T,as.is=T)
p.hwe<- bd$"P"[bd$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("HCGCON_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="HCGCON batch QQ-plot")
abline(a=0,b=1,lty=2)
bd.out <- subset(bd,bd$TEST=="ALL" & bd$P<=0.00001)$SNP
#hist(p.hwe,main="HWE test P value",xlab="P value")

be <-read.table(Fbe, header=T,as.is=T)
p.hwe<- be$"P"[be$"TEST"=="ALL"]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("HCPOPGEN_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="HCPOPGEN batch QQ-plot")
abline(a=0,b=1,lty=2)
be.out <- subset(be,be$TEST=="ALL" & be$P<=0.00001)$SNP
#hist(p.hwe,main="HWE test P value",xlab="P value")

bc <-read.table(Fbc, header=T,as.is=T)
p.hwe<- sort(bc$"P"[bc$"TEST"=="ALL"])
p.hwe<- p.hwe[p.hwe >= 1e-60]
x <- sort(-log10(ppoints(p.hwe))) #p.hwe pvalue list
y <- sort(-log10(p.hwe))
jpeg("HCDEMICK_hwe.jpg", res=300, width=10, height=5, units="in")
plot(x,y,ylab="Observed",xlab="Expected",main="HCDEMICK batch QQ-plot")
abline(a=0,b=1,lty=2)
bc.out <- subset(bc,bc$TEST=="ALL" & bc$P<= 0.00001)$SNP
#hist(p.hwe,main="HWE test P value",xlab="P value")

controlfails <-unique(c(ba.out,bb.out,bc.out,bd.out,be.out))
#cntrldup <- controlfails[duplicated(controlfails)]
fails <-c(case.out,ba.out,bb.out,bc.out,bd.out,be.out)
fails <- unique(fails)
write.table(fails, "SNPs_fail_HWE_5.txt", sep="\n",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(case.out, "SNPs_fail_case_HWE_5.txt", sep="\n",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(controlfails, "SNPs_fail_cntrls_HWE_5.txt", sep="\n",quote=FALSE,row.names=FALSE,col.names=FALSE)
#inter <-Reduce(intersect, list(case.out,controlfails))#=23
#cntrlinter <- Reduce(intersect, list(ba.out,bb.out,bc.out,bd.out,be.out))
cat(paste("Graphs and fail SNP list was successfully generated.\n",sep=""))

library(Vennerable)
v <- Venn(list(ba.out,bb.out,bc.out,bd.out,be.out,case.out),SetNames=c("HCDEBONN","HCDEKORA2","HCDEMICK","HCGCON","HCPOPGEN","CASES"))
jpeg("VennDiagram.jpg", res=300, width=10, height=5, units="in")
plot(v, doWeights = FALSE)#,type="AWFE") #ChowRuskey,AWFE , doEuler=TRUE)#
cat(paste("Venn diagram was successfully generated.\n",sep=""))
