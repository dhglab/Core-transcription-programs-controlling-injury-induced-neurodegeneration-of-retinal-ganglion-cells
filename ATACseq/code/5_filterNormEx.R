##***************************************************
# Process of peakCount matrix - two batches combined
# S1-8 and 9-20 combined
# 1. Preprocess Cov
# 2. filter Ex and peakAnno
# 3. cqn + deseq2 normalize 
# 4. remove batch effects 
# 5. remove outliers
# Yuyan 11-01-20
##****************************************************



library(DESeq2)
library(cqn)
library(limma)
library(ggplot2)
library(cowplot)
library(dplyr)
library(dummies)
library(corrplot)
library(gridExtra)
library(reshape2)
library(matrixStats)
library(pheatmap)
library(Hmisc)
library(RColorBrewer)

source("~/00_Utils.R")
dir='/u/home/y/ycheng41/nobackup-dhg/He_RGC/ATAC/CTCF'
indir=paste0(dir, '/nonTarget/DAR')
setwd(indir)
outdir=paste0(indir,"/0_Preprocess")
dir.create(outdir)

## 1. Preprocess Cov
AlignCov <- read.table("./AlignCOV.txt", sep = '\t',header = T,row.names = 1)
BioCov <- read.table(paste0(dir, "/BioCov.txt"),sep = '\t',header=T, row.names = 1)
BioCov <- BioCov[BioCov$Injection == 'nonTarget',]

BioCov <- BioCov[rownames(AlignCov),]

Cov <- cbind(AlignCov,BioCov)
Cov$Time <- factor(Cov$Time, levels = c(0,1,3))
Cov$Batch <- as.factor(Cov$Batch)
## 2. filter Ex and peakAnno
EXPR <- read.table("./peaks_countMatrix.txt", header = TRUE, row.names = 1)
Ex <- EXPR[,-c(1:5)]
colnames(Ex) <- gsub('X.u.home.y.ycheng41.nobackup.dhg.He_RGC.ATAC.CTCF.nonTarget.Bowtie2.','',colnames(Ex))
colnames(Ex) <- gsub('.unique.sorted.rmDup.bam','',colnames(Ex))

Ex <- Ex[,rownames(Cov)]
peakAnno <- read.table("./peak.gc.txt",sep = '\t')
colnames(peakAnno) <- c("Chr","Start","End", "pctGC", "Length")

peakAnno$Start <- peakAnno$Start + 1 #featureCounts is 1-based

#all((peakAnno$Start) == EXPR$Start); all(peakAnno$End == EXPR$End)
rownames(peakAnno) <- rownames(Ex)

# filter Ex for low-count peaks
idx <- rowSums((Ex)>5) >= ncol(Ex)
Ex <- Ex[idx,]
peakAnno <- peakAnno[idx,]

save(peakAnno,Ex,Cov,file=paste0(outdir, "/filterEx.rda"))

##3. cqn + deseq2 normalize 


load(paste0(indir,'/0_Preprocess/filterEx.rda'))

# run cqn to generate the normalization factors matrix to account for GC and length bias
cqnObject=cqn(counts = Ex,x=peakAnno$pctGC, lengths = peakAnno$Length, 
              sizeFactors =Cov$UniqMappedRead)
cqnOffset <- cqnObject$glm.offset
cqnNormFactors <- exp(cqnOffset)
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))


#deseq2 normalization
dds <- DESeqDataSetFromMatrix(countData = Ex, colData = Cov, design = ~ 1)
normalizationFactors(dds) <- normFactors#dds <- estimateSizeFactors(dds)
normEx <- counts(dds, normalized=T)

normEx <- log2(normEx) # have compared VST vs log2 - similar
save(normEx,peakAnno,Cov, file = paste0(outdir, '/normEx.rda'))

#Corrplot of pca with covariates
Cov.sub <- Cov[,c("DuplicateRate","UniqMappedRead","FRiPNarrowPeak","Batch")]

Cov.sub$Batch <- as.numeric(Cov.sub$Batch)
Cov1 <- Cov.sub[,!unlist(lapply(Cov.sub,is.factor))] #numeric covariates
Cov2 <- Cov.sub[,unlist(lapply(Cov.sub,is.factor)),drop=F] 
Dummy <- NULL #use dummy table to change factoric covariates
if (dim(Cov2)[2] >0) { Dummy <- dummy.data.frame(Cov2, names = colnames(Cov2)); COVa <- cbind(Cov1, Dummy)}else{
  COVa <- Cov1
}




#PCA on normEx
num=c(1000, 5000,nrow(normEx))[3]
rv = rowVars(as.matrix(normEx)) ;  select = order(rv, decreasing=TRUE)[seq_len(num)]

pca = prcomp(t(normEx[select,]))
summary(pca)
pcaData=as.data.frame(pca$x)

prcomp_subset <- pcaData[,c(1:5)]


COVB <- cbind(prcomp_subset, COVa)
COVB <- COVB[sapply(COVB,is.numeric)]

res <- rcorr(as.matrix(COVB)) #Hmisc methods to calculate both correlation


pdf(file=paste0(outdir,'/', "normExPCA_corrplot.pdf"),width = 8,height = 8)
par(mar=c(2,2,1,1))
corrplot(res$r, p.mat=res$P, sig.level = 0.005, insig="label_sig", diag = FALSE, type = 'lower')#,upper="ellipse",title="Corrplot")
dev.off()


## 4. remove batch effects, as it strongly correlates with PC1
regEx <- removeBatchEffect(normEx, Cov$Batch)
save(regEx, Cov, peakAnno, file=paste0(outdir,'/regressEx.rda'))

#PCA on regEx
num=c(1000, 5000,nrow(regEx))[3]
rv = rowVars(as.matrix(regEx)) ;  select = order(rv, decreasing=TRUE)[seq_len(num)]

pca = prcomp(t(regEx[select,]))
summary(pca)
pcaData=as.data.frame(pca$x)

prcomp_subset <- pcaData[,c(1:5)]
#pcaResult.cqn<-prcomp(t(RPKM.cqn)) ; prcomp_subset.cqn <- pcaResult.cqn$x[,c(1:5)] ; 
#colnames(prcomp_subset.cqn) <- gsub("PC",paste0("cqn.PC"),colnames(prcomp_subset.cqn))


COVB <- cbind(prcomp_subset, COVa)
COVB <- COVB[sapply(COVB,is.numeric)]

res <- rcorr(as.matrix(COVB)) #Hmisc methods to calculate both correlation


pdf(file=paste0(outdir,'/',"regExPCA_corrplot.pdf"),width = 8,height = 8)
par(mar=c(2,2,1,1))
corrplot(res$r, p.mat=res$P, insig="label_sig", diag = FALSE, type = 'lower')#,upper="ellipse",title="Corrplot")
dev.off()


## 5. sample distance
Cov$Group <- paste0(Cov$Injection,'_day',Cov$Time)
Cov$Group <- factor(Cov$Group, levels = c("nonTarget_day0","nonTarget_day1","nonTarget_day3",
                                          "CTCF_day0","CTCF_day1", "CTCF_day3"))
annotation <- data.frame(row.names = rownames(Cov),Cov$Injection,Cov$Treatment)

colnames(annotation) <- c('Injection','Treatment')
sampleDists <- dist(t(regEx))
plotEx <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file=paste0(outdir, '/reg.sampleCor.pdf'), width = 6,height = 6)
par(mar=c(5,3,3,5))
pheatmap(plotEx,
         show_rownames = T,
         show_colnames = T,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         annotation_col = annotation,silent = F)
dev.off()





