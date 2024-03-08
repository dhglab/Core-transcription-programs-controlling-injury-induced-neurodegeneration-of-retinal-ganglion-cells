## load R pacakges ##
liblist <- c("edgeR", "readxl", "dplyr","fastDummies","Hmisc","corrplot",
             "matrixStats", "ggplot2","RColorBrewer","cowplot", "ComplexHeatmap",
             "circlize",'ggpubr')
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = FALSE)))

            
## util functions ##

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

cal_cov <- function(AlignCov,outdir){
  
  
  #PCA on sequencing PC
  
  pca = prcomp(AlignCov)
  cut <- 0.01
  index <- sum(summary(pca)$importance[2,] > cut)# how many PC
  index2 <- max(2,index)
  pdf(paste0(outdir, '/AlignCov_scree_PCA1.pdf'), h=6,w=8)
  barplot(summary(pca)$importance[2,],xlab="Principle Component",ylab="Variance Explained")
  abline(h=cut,col="blue")
  dev.off()
  
  seqPC <- pca$x[,1:index2]
  colnames(seqPC) <- gsub("PC","seqPC",colnames(seqPC))
  
  
  return(seqPC)
}  
cal_cov_ex_corr <- function(Ex,COV){
  num=c(1000, 5000,nrow(Ex))[3]
  rv = rowVars(as.matrix(Ex)) ;  select = order(rv, decreasing=TRUE)[seq_len(num)]
  
  pca = prcomp(t(Ex[select,]))
  summary(pca)
  pcaData=as.data.frame(pca$x)
  
  prcomp_subset <- pcaData[,c(1:5)]
  
  
  COVB <- cbind(prcomp_subset,COV)
  COVB <- COVB[sapply(COVB,is.numeric)]
  
  res <- rcorr(as.matrix(COVB)) #Hmisc methods to calculate both correlation
  return(res)
}

WGCNA.outlier <- function(Ex, outdir, name_prefix){
  require(WGCNA); require(ggplot2)
  normalized_adj <- (0.5 + 0.5*bicor(t(Ex)))^2  ## normalizing expression levels, because the function only takes value between 0-1
  netsummary <- fundamentalNetworkConcepts(normalized_adj)
  connectivity <- netsummary$Connectivity
  connectivity.zscore <- (connectivity-mean(connectivity))/sqrt(var(connectivity))
  connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
  
  pdf(paste0(outdir, "/", name_prefix, "_sample_connectivity.pdf"), height = 10, width = 10)
  p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
  p <- p+ geom_hline(aes(yintercept = -2))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
  print(p)
  dev.off()
  
}
lm.fun <- function(Expr,expression,designvars) {
  # but study should be in the env - as should "Method" and "Add"
  # expr = the formula, design vars is a dataframe of
  library(nlme) ; library(Hmisc) ; library(lme4)
  #print(info(Cov)) ;  print(info(Expr)) ; print(info(designvars)) ; print(expression)
  
  residuals1 <- residuals <- matrix(NA,nrow=nrow(Expr),ncol=ncol(Expr))
  gene.assoc.post <- NULL ; gene.assoc <- NULL
  
  for (i in 1:nrow(Expr)) {
    # notes, chris loged age and rin
    
    if (i %% 500 == 0) {cat(print(paste("On gene ",i,"\n",sep="")))}
    thisExpr <- as.numeric(Expr[i,])
    designmatrix <- data.frame(thisExpr, designvars)
    m1 <- eval(parse(text=expression));   # adding the +1 really changes the resilts
    # get residuals + intercept and tissue
    residuals[i,] <-  summary(m1)$coefficients[1] + residuals(m1)
    if ( is.null(gene.assoc) ) {
      gene.assoc <- matrix(0, nrow=nrow(Expr), ncol=nrow(summary(m1)$coefficients))
      gene.assoc.post <- 0 * gene.assoc  }
    gene.assoc[i,] <- summary(m1)$coefficients[,3]^2  # Chi-squared stat
    # this has the association of each gene to each covariate
    
    # get post-regression covariates - what does this mean, association of gene to covariate AFTER regression - so doing lmr on the newly regressed dataset
    df <- designvars
    df$thisExpr <- residuals[i,]
    m1.post <- eval(parse(text=paste0(strsplit(expression,", data")[[1]][1],", data=df)")))
    gene.assoc.post[i,] <- summary(m1.post)$coefficients[,3]^2  # Chi-squared stat
    ###### the poast assoc scores are similar for tissues, BUT alot smaller(or significant) for the other factors. emplies regressed out (or more significant)
  }
  colnames(residuals) <- colnames(Expr)
  rownames(residuals) <- rownames(gene.assoc) <- rownames(gene.assoc.post) <- rownames(Expr)
  colnames(gene.assoc) <- colnames(gene.assoc.post) <- rownames(summary(m1.post)$coefficients)
  print(info(residuals))   ; print(info(Expr))
  rownames(residuals) <- rownames(gene.assoc) <- rownames(gene.assoc.post) <- rownames(Expr)
  colnames(gene.assoc) <- colnames(gene.assoc.post) <- rownames(summary(m1.post)$coefficients)
  print(info(residuals))   ; print(info(Expr))
  
  return(list(as.data.frame(residuals), gene.assoc, gene.assoc.post))
}
edgeR.process <- function(Cov, Ex, geneAnno){
  require(edgeR)
  y <- DGEList(counts=Ex, samples=Cov, genes=geneAnno)
  keep <- filterByExpr(y,min.count = 5, min.prop = 0.5)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  #normalize by tmm
  y <- calcNormFactors(y)
  
  return(y)
  
}

edgeR.dex <- function(this.tf, y, design, outdir){
  require(edgeR)
  mod <- model.matrix(as.formula(design), data=y$samples)
  colnames(mod) <- gsub('group','',colnames(mod))
  
  y <- estimateDisp(y, mod, robust=TRUE)
  
  pdf(paste0(outdir, '/', this.tf,'.dispersion.pdf'), width = 5,height = 5)
  plotBCV(y)
  dev.off()
  
  fit <- glmQLFit(y, mod)
  
  DEsets <- fitsets <- list()
  cmd1 <- paste0("con1 <- makeContrasts(", this.tf,'.0d_vs_non.0d=', this.tf, '.0d-non.0d,', 'levels = mod)')
  eval(parse(text = cmd1))
  cmd2 <- paste0("con2 <- makeContrasts(", this.tf,'.3d_vs_non.3d=', this.tf, '.3d-non.3d,', 'levels = mod)')
  eval(parse(text = cmd2))
  
  
  #con1 <- paste0(this.tf,'.0d_vs_non.0d=', this.tf, '.0d-non.0d')
  #con2 <- paste0(this.tf,'.3d_vs_non.3d=', this.tf, '.3d-non.3d')
  
  #fitsets[[1]] <- glmQLFTest(fit, contrast=makeContrasts(con1, levels=mod))
  # fitsets[[2]] <- glmQLFTest(fit, contrast=makeContrasts(con2, levels=mod))
  
  fitsets[[1]] <- glmQLFTest(fit, contrast=con1)
  fitsets[[2]] <- glmQLFTest(fit, contrast=con2)
  
  
  DEsets <- lapply(fitsets, function(i){topTags(i,n = Inf, adjust.method = 'fdr',sort.by = 'none')})
  names(DEsets) <- names(fitsets) <- c(paste0(this.tf,'.0d_vs_non.0d'),paste0(this.tf,'.3d_vs_non.3d'))
  
  save(DEsets, fitsets,file=paste0(outdir,'/', this.tf, '.edgeR.rda'))
}

NumDEX <- function(fitsets, outdir, prefix){
  decide <- matrix(c("fdr",0.05, "fdr", 0.1, "none", 0.005, "none", 0.01), nrow = 4, ncol=2, byrow = TRUE)
  mysum <- as.list(1:nrow(decide))
  #mynum <- 0
  maxmax <- 0
  index<-0
  for (test in 1:nrow(decide)){
    for(fit in fitsets){
      index <- index +1
      results <- decideTests(fit, adjust.method = decide[test,1],
                             p = as.numeric(decide[test,2]))
      mysum[[index]] <- summary(results)[,1] # this is a list of number for -1, 0, 1
      maxmax <- max(c(maxmax, as.vector(mysum[[index]][c(1,3)])))
    }
  }
  
  fitset_names <- names(fitsets)
  plotMe1<-numeric(length = length(fitset_names))
  plotMe2<-numeric(length = length(fitset_names))
  
  
  pdf(paste0(outdir,'/', prefix, ".NumDEX.pdf"),width = 8, height = 2.5)
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
}
