rm(list=ls())


source("./code/00_utils.R")

tfs <- c('ATF3','ATF4','CHOP','CEBPG')
load('./output/preprocess.rda')
outdir <- './output'

## identify sample outliers
# WGCNA sample connectivity absolute z score > 2

WGCNA.outlier(t(Ex), outdir, 'raw') # identified ATF4.3d.5, ATF4.3d.6

## remove ATF4.3d.5, ATF4.3d.6
toDel <- c('ATF4.3d.5','ATF4.3d.6')
idx <- is.na(match(rownames(Cov), toDel))
Cov1 <- Cov[idx,]
Ex1 <- Ex[, rownames(Cov1)]

## normalize by edgeR TMM
y <- edgeR.process(Cov1, Ex1, geneAnno)
normEx <- cpm(y, log=F)

##  regress seqPC1, seqPC2, seqPC3, seqPC4

ToKeep = c("seqPC1", "seqPC2", "seqPC3", "seqPC4")
designvars <- Cov1[,ToKeep,drop=F]
expression <- "lm(thisExpr ~ seqPC1 + seqPC2 + seqPC3 + seqPC4, data = designmatrix)"
Method <- "lm"

LIST <- lm.fun(normEx,expression,designvars)
regEx <- LIST[[1]]

#check outlier from normzlied, seqPC regressed Ex
WGCNA.outlier(t(regEx), outdir, 'rm1_reg') # identified CEBPG.3d.7

### remove ATF4.3d.5, ATF4.3d.6, 'CEBPG.3d.7'

toDel <- c('ATF4.3d.5','ATF4.3d.6','CEBPG.3d.7')
idx <- is.na(match(rownames(Cov), toDel))
Cov1 <- Cov[idx,]
Ex1 <- Ex[, rownames(Cov1)]

y <- edgeR.process(Cov1, Ex1, geneAnno)
normEx <- cpm(y, log=T)


WGCNA.outlier(t(normEx), outdir, 'rm2_norm')

ToKeep = c("seqPC1", "seqPC2", "seqPC3", "seqPC4")
designvars <- Cov1[,ToKeep,drop=F]
expression <- "lm(thisExpr ~ seqPC1 + seqPC2 + seqPC3 + seqPC4, data = designmatrix)"
Method <- "lm"
LIST <- lm.fun(normEx,expression,designvars)
regEx <- LIST[[1]]
WGCNA.outlier(t(regEx), outdir, 'rm2_reg')


# sample distance heatmap - further removing "non.3d.9","non.3d.11", 'ATF4.3d.7'


  annotation <- data.frame(row.names = rownames(Cov),Cov$Target)
  
  colnames(annotation) <- c('Treatment')
  sampleDists <- dist(t(normEx))
  plotEx <- as.matrix( sampleDists )
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(file=paste0(dir, '/normEx_sampleDist.pdf'), width = 6,height = 6)
  par(mar=c(5,3,3,5))
  pheatmap(plotEx,
           show_rownames = F,
           show_colnames = T,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           annotation_col = annotation,silent = F)
  dev.off()
  

