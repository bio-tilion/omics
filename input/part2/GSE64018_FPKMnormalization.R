rm(list=ls())
options(stringsAsFactors=FALSE)

## Depending on how the files are loaded, some of the code may need to be changed
datExpr <- read.table("../datCounts/countlevel_12asd_12ctl.txt",sep="\t",header=TRUE)
datMeta <- read.table("../datMeta/12v12_GEO_upload.txt",sep="\t",header=TRUE)

datFPM <- sweep(log2(datExpr + 1), 2, log2(apply(datExpr,2,sum)/10^6))
keep <- apply(datFPM>1,1,sum)>(0.80*ncol(datFPM)) ## Keep those with > 1 FPM in 80% of samples
datExpr <- datFPM[keep,]

## Do a PCA of the sequencing statistics of the full sample
colnames(datMeta) <- gsub("characteristics..","",colnames(datMeta))
datSeq <- datMeta[,c(16:30)] ## All 15 columns, "TotalReads.picard" to "PropExonicReads.HTSC"
datSeqNorm <- t(scale(datSeq,scale=F))
PC.datSeq <- prcomp(datSeqNorm);
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
print(varexp[1:2])
topPC.datSeq <- PC.datSeq$rotation[,1:2]; ## Explains 99% of variance in datSeq
colnames(topPC.datSeq) <- c("SeqPC1 - Depth","SeqPC2 - GC/Length")

## Get the data
condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
age <- as.numeric(datMeta[,"Age"])
sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
region <- as.numeric(as.factor(datMeta[,"RegionID"]))-1
RIN <- as.numeric(datMeta[,"RIN"])
bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
seqStatPC1 <- topPC.datSeq[,1]
seqStatPC2 <- topPC.datSeq[,2]

varnames <- c("condition","age","sex","RIN","bank","seqStatPC1","seqStatPC2")
Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datExpr),ncol=7)
colnames(Bmat) <- paste("beta",varnames,sep=".")
colnames(SEmat) <- paste("SE",varnames,sep=".")
colnames(Pmat) <- paste("p",varnames,sep=".")
rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datExpr)

## Adjusted values
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)
regvars <- data.frame(condition=condition,age=age,sex=sex,RIN=RIN,bank=bank,seqStatPC1=seqStatPC1,seqStatPC2=seqStatPC2)
coefmat <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(regvars)+1)

for (i in 1:nrow(datExpr)) {
  if (i %% 1000 == 0) {cat(paste("On gene ",i,"\n",sep=""))}
  thisExpr <- as.numeric(datExpr[i,])
  lm.out <- lm(thisExpr ~ condition + age + sex + RIN + bank + seqStatPC1 + seqStatPC2)

  coef <- coef(lm.out) ## Get the coefficients from the model
  coefmat[i,] <- coef
  datExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lm.out$residuals 
  
  tabOut <- summary(lm.out)$coefficients
  Bmat[i,] <- tabOut[-c(1),"Estimate"]
  SEmat[i,] <- tabOut[-c(1),"Std. Error"]
  Pmat[i,] <- tabOut[-c(1),"Pr(>|t|)"]
}

## Write adjusted expression values
write.table(datExpr.reg,"../datAdjusted/adjfpkm_12asd_12ctl.txt",sep="\t")

## Get the information for SRRM4
pdf("./SRRM4_expression.pdf",height=8,width=8)
par(mfrow=c(2,2))
tResults <- t.test(as.numeric(datExpr["ENSG00000139767",])~as.factor(datMeta[,"Diagnosis"]))
boxplot(as.numeric(datExpr["ENSG00000139767",])~as.factor(datMeta[,"Diagnosis"]),ylab="Normalized FPKM",main="SRRM4 - unadjusted expression",cex.main=1,cex.lab=1,cex.axis=1)
points(y=as.numeric(datExpr["ENSG00000139767",]),x=as.numeric(as.factor(datMeta[,"Diagnosis"])),pch=19,col="orange")

tResults <- t.test(datExpr.reg["ENSG00000139767",]~as.factor(datMeta[,"Diagnosis"]))
boxplot(datExpr.reg["ENSG00000139767",]~as.factor(datMeta[,"Diagnosis"]),ylab="Normalized FPKM",main="SRRM4 - adjusted expression",cex.main=1,cex.lab=1,cex.axis=1)
points(y=datExpr.reg["ENSG00000139767",],x=as.numeric(as.factor(datMeta[,"Diagnosis"])),pch=19,col="orange")

dev.off()
