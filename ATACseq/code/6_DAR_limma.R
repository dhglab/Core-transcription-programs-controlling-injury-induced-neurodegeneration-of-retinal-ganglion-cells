##***************************************************
# differential expression on nonTarget Ex S1-26
# Yuyan 11-01-20
# module load R/4.0.2
##****************************************************
library(limma)
library(ggpubr)
library(dplyr)
library(cowplot)
library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)
library(Vennerable)
source("~/00_Utils_common.R")

dir='/u/home/y/ycheng41/nobackup-dhg/He_RGC/ATAC/CTCF'
indir=paste0(dir, '/nonTarget/DAR')
setwd(indir)

outdir=paste0(indir,"/1_limma")
dir.create(outdir)
load(paste0(dir,'/DAR/0_Preprocess/normEx.rda'))


Cov$Batch <- as.factor(Cov$Batch)
design <- (~Cov$Time + Cov$Batch)
mod <- model.matrix(as.formula(design))
colnames(mod) <- gsub('Cov\\$','',colnames(mod))

## 1. differential analysis

# a. limma
fit <- lmFit(normEx,mod)
colnames(mod)
#"(Intercept)"    "Time1"          "Time3"          "FRiPNarrowPeak"
##intercept: Time0
##Time1: Time1 - Time0
##Time3: Time3 - Time0
DEsets <- fitsets <- list()

fitsets[[1]] <- eBayes(contrasts.fit(fit, contrasts = c(0,1,0,0)))
fitsets[[2]] <- eBayes(contrasts.fit(fit, contrasts = c(0,0,1,0)))
fitsets[[3]] <- eBayes(contrasts.fit(fit, contrasts = c(0,1,1,0)))

DEsets <- lapply(fitsets, function(i){topTable(i,number = Inf, adjust.method = 'fdr',sort.by = 'none')})
names(DEsets) <- names(fitsets) <- c("nonTarget_Day1vsDay0","nonTarget_Day3vsDay0","nonTarget_Day3vsDay1")

DAR <- do.call(cbind, DEsets)
DAR.logFC <- DAR[,grep("logFC",colnames(DAR))]
DAR.p <- DAR[,grep(paste(c("P.Value", "adj.P.Val"),collapse = "|"),colnames(DAR))]


DAR.AveExpr <- DAR[,grep("AveExpr",colnames(DAR))]
DAR <- cbind(DAR.logFC,DAR.p,DAR.AveExpr)

save(DEsets, fitsets,file=paste0(outdir,'/limma.rda'))



#b.number of DARs 

decide <- matrix(c("fdr",0.05, "fdr", 0.1, "none", 0.005, "none", 0.01), nrow = 4, ncol=2, byrow = TRUE)
mysum <- as.list(1:nrow(decide))
#mynum <- 0
maxmax <- 0
index<-0
for (test in 1:nrow(decide)){
  for(fit in fitsets[1:2]){
    index <- index +1
    results <- decideTests(fit, adjust.method = decide[test,1],
                           p = as.numeric(decide[test,2]))
    mysum[[index]] <- summary(results)[,1] # this is a list of number for -1, 0, 1
    maxmax <- max(c(maxmax, as.vector(mysum[[index]][c(1,3)])))
  }
}

fitset_names <- names(DEsets)[1:2]
plotMe1<-numeric(length = length(fitset_names))
plotMe2<-numeric(length = length(fitset_names))


pdf(paste0(outdir, "/NumDAR.pdf"),width = 8, height = 2.5)
par(mfrow=c(1,nrow(decide)))


index <- 0
yy <- NULL
for(test in 1:nrow(decide)){
  for(i in 1:length(fitset_names)){
    index <- index + 1
    plotMe1[i] <- as.numeric(as.vector(mysum[[index]][3]))
    plotMe2[i] <- as.numeric(as.vector(mysum[[index]][1]))
  }
  
  maxData1 = max(plotMe1)
  maxData2 = max(plotMe2)
  yy <- barplot(plotMe1, horiz = TRUE, col = "#B31B21", xlim = c(-maxmax*1.2, maxmax*1.2), main = paste("Gene Changes \np<", decide[test,2], ",",
                                                                                                    decide[test,1], sep =""),cex.main=0.9)
  yy <- barplot(-plotMe2, horiz = TRUE, col = "#1465AC", add = T, axes = FALSE)
  #xx<-vector("integer",length(fitset_names))
  text(yy,fitset_names,cex=0.9)
  text((plotMe1+10)*0 + 0.9*maxData1,yy+0.2,format(plotMe1,digits=3),cex = 0.9)
  text((-plotMe2-10)*0 - 0.9*maxData2,yy+0.2,format(plotMe2,digits=3), cex = 0.9)
}

dev.off()



#c. topDAR
topDAR = data.frame(matrix(NA, 0, ncol(DAR), dimnames = list(c(), colnames(DAR))))

for(i in 1:nrow(DAR)){
  for(j in seq(2,ncol(DAR.p),2)){
    if(!is.na(DAR.p[i, j])){
      if(DAR.p[i, j] < 0.1){
        topDAR = rbind(topDAR, DAR[i,])
        break
      }
    }
  }
}

save(DAR, topDAR, peakAnno, Cov, file=paste0(outdir,'/DAR.rda'))


#d. MA plot - use customized ggmaplot2.R
source("~/apps/ggpubr/00_Utils.R")
source("~/apps/ggpubr/ggmaplot2.R")

myplot1 <-  DEsets[[1]][,c("AveExpr",'logFC','adj.P.Val')]
colnames(myplot1) <- c('baseMeanLog2', 'log2FoldChange', 'padj')
up.n1=100*table(myplot1$log2FoldChange > 0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
down.n1=100*table(myplot1$log2FoldChange < -0.3 &  myplot1$padj < 0.1)[2]/nrow(myplot1)
up.n1=as.numeric(format(round(up.n1, 2), nsmall = 2))
down.n1=as.numeric(format(round(down.n1, 2), nsmall = 2))

p1 <- ggmaplot2(myplot1, main = paste0("Day 1 vs Day 0\nUp: ", up.n1, "%; Down: ", down.n1, '%'),
         xlab = 'Log2 mean accessibility', fdr = 0.1, fc = 1.3, size = 0.005,alpha=0.4,
         palette = c("darkgrey", "#1465AC","#B31B21"),
         genenames = rep('',nrow(myplot1)),
         legend = "none", 
         font.y= c("plain", 9),
         font.x= c("plain", 9),
         font.main = c("plain", 9),
         font.tickslab = c(8, "plain"),
         ylim=c(-2.5,2.5),xlim=c(0,16),
         ggtheme = ggplot2::theme_minimal())


myplot2 <-  DEsets[[2]][,c("AveExpr",'logFC','adj.P.Val')]
colnames(myplot2) <- c('baseMeanLog2', 'log2FoldChange', 'padj')
up.n2=100*table(myplot2$log2FoldChange > 0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
up.n2=as.numeric(format(round(up.n2, 2), nsmall = 2))
down.n2=100*table(myplot2$log2FoldChange < -0.3 &  myplot2$padj < 0.1)[2]/nrow(myplot2)
down.n2=as.numeric(format(round(down.n2, 2), nsmall = 2))

p2 <- ggmaplot2(myplot2,main = paste0("Day 3 vs Day 0\nUp: ", up.n2, "%; Down: ", down.n2, '%'),
                xlab = 'Log2 mean accessibility', fdr = 0.1, fc = 1.3, size = 0.005,alpha=0.4,
               palette = c("darkgrey", "#1465AC","#B31B21"),
               genenames = rep('',nrow(myplot2)),
               legend = "none",
               font.x= c("plain", 9),
               font.y= c("plain", 9),
               font.main = c("plain", 9),
               font.tickslab = c(8, "plain"),
               ylim=c(-2.5,2.5),xlim=c(0,16),
               ggtheme = ggplot2::theme_minimal())
png(file=paste0(outdir,'/MAplotDAR.png'), width = 4.5,height = 2, units = 'in',res = 600)
plot_grid(p1,p2,  nrow = 1,labels = NULL)
dev.off()

