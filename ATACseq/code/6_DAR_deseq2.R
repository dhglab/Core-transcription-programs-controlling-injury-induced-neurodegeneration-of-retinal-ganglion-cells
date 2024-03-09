library(DESeq2)
library(cqn)
library(cowplot)
library(ggpubr)
source("~/apps/ggpubr/00_Utils.R")
source("~/apps/ggpubr/ggmaplot2.R")


dir='/u/home/y/ycheng41/nobackup-dhg/He_RGC/ATAC/CTCF'
indir=paste0(dir, '/nonTarget/DAR')
outdir=paste0(indir,"/1_limma")


load(paste0(indir,'/0_Preprocess/filterEx.rda'))

# run cqn to generate the normalization factors matrix to account for GC and length bias
cqnObject=cqn(counts = Ex,x=peakAnno$pctGC, lengths = peakAnno$Length, 
              sizeFactors =Cov$UniqMappedRead)
cqnOffset <- cqnObject$glm.offset
cqnNormFactors <- exp(cqnOffset)
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))


dds <- DESeqDataSetFromMatrix(countData = Ex, colData = Cov, design = ~ 1)
design(dds) <-  ~Time + Batch

normalizationFactors(dds) <- normFactors
dds <- DESeq(dds)

DEsets <- list()

DEsets[[1]] <- results(dds, contrast=c('Time',1,0), pAdjustMethod='fdr')
DEsets[[2]] <- results(dds, contrast=c('Time',3,0), pAdjustMethod='fdr')
DEsets[[3]] <- results(dds, contrast=c('Time',3,1), pAdjustMethod='fdr')

names(DEsets) <-  c("nonTarget_Day1vsDay0","nonTarget_Day3vsDay0","nonTarget_Day3vsDay1")

DAR <- do.call(cbind, DEsets)
DAR.logFC <- DAR[,grep("log2FoldChange",colnames(DAR))]
DAR.p <- DAR[,grep(paste(c("pvalue", "padj"),collapse = "|"),colnames(DAR))]
DAR <- cbind(DAR.logFC, DAR.p)

myplot1 <-  DEsets[[1]][,c("baseMean",'log2FoldChange','padj')]
#up.n1=100*table(myplot1$log2FoldChange > 0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
#down.n1=100*table(myplot1$log2FoldChange < -0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
up.n1=4.11#as.numeric(format(round(up.n1, 2), nsmall = 2))
down.n1=5.14#as.numeric(format(round(down.n1, 2), nsmall = 2))

p1 <- ggmaplot2(myplot1, main = paste0("Day1 vs Day0\nUp: ", up.n1, "%; Down: ", down.n1, '%'),
                xlab = 'Log2 mean accessbility', fdr = 0.1, fc = 1.3, size = 0.0005,alpha=0.4,xlim=c(0,16),ylim=c(-3,3),
                palette = c("darkgrey", "#1465AC","#B31B21"),
                genenames = rep('',nrow(myplot1)),
                legend = "none", 
                font.y= c("plain", 9),
                font.x= c("plain", 9),
                font.main = c("plain", 9),
                font.tickslab = c(8, "plain"),
                ggtheme = ggplot2::theme_minimal())


myplot2 <-  DEsets[[2]][,c("baseMean",'log2FoldChange','padj')]
#up.n2=100*table(myplot2$log2FoldChange > 0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
up.n2=14.65#as.numeric(format(round(up.n2, 2), nsmall = 2))
#down.n2=100*table(myplot2$log2FoldChange < -0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
down.n2=12.6#as.numeric(format(round(down.n2, 2), nsmall = 2))

p2 <- ggmaplot2(myplot2,main = paste0("Day3 vs Day0\nUp: ", up.n2, "%; Down: ", down.n2, '%'),
                xlab = 'Log2 mean accessbility', fdr = 0.1, fc = 1.3, size = 0.0005,alpha=0.4,xlim=c(0,16),ylim=c(-3,3),
                palette = c("darkgrey", "#1465AC","#B31B21"),
                genenames = rep('',nrow(myplot2)),
                legend = "none",
                font.x= c("plain", 9),
                font.y= c("plain", 9),
                font.main = c("plain", 9),
                font.tickslab = c(8, "plain"),
                ggtheme = ggplot2::theme_minimal())


png(file=paste0(outdir,'/MAplotDAR.png'), width = 4.5,height = 2, units = 'in',res = 600)
plot_grid(p1,p2,  nrow = 1,labels = NULL)
dev.off()


load(paste0(dir,'/nonTarget/Ex_all.rda'))
DEX <- read.table(paste0(dir, "/nonTarget/DEX.txt"), sep = '\t',header = T)

ave.expr1 <- rnaEx[,rnaCov$Time==0 | rnaCov$Time==1]
ave.expr1 <- rowMeans(ave.expr1)
ave.expr2 <- rnaEx[,rnaCov$Time==0 | rnaCov$Time==3]
ave.expr2 <- rowMeans(ave.expr2)


myplot1 <-  cbind(ave.expr1, DEX[,c(1,5)])
colnames(myplot1) <- c('baseMeanLog2', 'log2FoldChange', 'padj')

up.n1=100*table(myplot1$log2FoldChange > 0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
down.n1=100*table(myplot1$log2FoldChange < -0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
up.n1=as.numeric(format(round(up.n1, 2), nsmall = 2))
down.n1=as.numeric(format(round(down.n1, 2), nsmall = 2))

p1 <- ggmaplot2(myplot1, main = paste0("Day1 vs Day0\nUp: ", up.n1, "%; Down: ", down.n1, '%'),
                fdr = 0.1, fc = 1.3, size = 0.01,alpha=0.7,
                palette = c("darkgrey", "#1465AC","#B31B21"),
                genenames = rep('',nrow(myplot1)),
                legend = "none", 
                font.y= c("plain", 9),
                font.x= c("plain", 9),
                font.main = c("plain", 9),
                font.tickslab = c(8, "plain"),
                ylim=c(-4,4),xlim=c(0,16),
                ggtheme = ggplot2::theme_minimal())


myplot2 <-  cbind(ave.expr2, DEX[,c(2,6)])
colnames(myplot2) <- c('baseMeanLog2', 'log2FoldChange', 'padj')

up.n2=100*table(myplot2$log2FoldChange > 0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
up.n2=as.numeric(format(round(up.n2, 2), nsmall = 2))
down.n2=100*table(myplot2$log2FoldChange < -0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
down.n2=as.numeric(format(round(down.n2, 2), nsmall = 2))

p2 <- ggmaplot2(myplot2, main = paste0("Day3 vs Day0\nUp: ", up.n2, "%; Down: ", down.n2, '%'),
                fdr = 0.1, fc = 1.3, size = 0.01,alpha=0.7,
                palette = c("darkgrey", "#1465AC","#B31B21"),
                genenames = rep('',nrow(myplot2)),
                legend = "none",
                font.x= c("plain", 9),
                font.y= c("plain", 9),
                font.main = c("plain", 9),
                font.tickslab = c(8, "plain"),
                ylim=c(-4,4),xlim=c(0,16),
                ggtheme = ggplot2::theme_minimal())
png(file=paste0(outdir,'/MAplotDEX.png'), width = 4.5,height = 2, units = 'in',res = 300)
plot_grid(p1,p2,  nrow = 1,labels = NULL)
dev.off()


