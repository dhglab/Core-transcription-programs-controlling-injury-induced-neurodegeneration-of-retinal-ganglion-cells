rm(list=ls())


source("./code/00_utils.R")

## identify covariates for regression

tfs <- c('ATF3','ATF4','CHOP','CEBPG')
load('../rawEx.rda')
outdir <- '../output'; if(!file.exists(outdir)){dir.create(outdir)}

# biological covariant
BioCov <- Cov[,c('TF','Batch','Sex'),drop=F]
BioCov <- lapply(BioCov, as.factor) %>% as.data.frame()
BioCov$Batch <- as.numeric(BioCov$Batch); BioCov$Sex <- as.numeric(BioCov$Sex)
Cov1 <- BioCov[,!unlist(lapply(BioCov,is.factor)), drop=F] #numeric covariates
Cov2 <- BioCov[,unlist(lapply(BioCov,is.factor)),drop=F] #factorial covariates
Dummy <- NULL #use dummy table to change factoric covariates
if (dim(Cov2)[2] >0) { Dummy <- dummy_cols(Cov2); COVa <- cbind(Cov1, Dummy)}else{
  COVa <- Cov1
}

#sequencing covariant
AlignCov <- Cov[,c(10:ncol(Cov))]
AlignCov <- t(AlignCov) %>% scale_rows() %>% t()

#PCA of sequencing covariates
COVb <- cal_cov(AlignCov,outdir)

#correlation between expression PCs, sequencing and biological covariates

AlignBioCov <- cbind(COVa, COVb)
res <- cal_cov_ex_corr(Ex, AlignBioCov) 

#identified seqPC1-4 to regress
#Batch and Sex are correlated with seqPCs

AlignBioCov <- AlignBioCov %>% select(Batch, Sex, seqPC1, seqPC2, seqPC3, seqPC4)
res <- cal_cov_ex_corr(Ex, AlignBioCov) 

pdf(file=paste0(outdir,'/cov_Ex_corrplot.pdf'),width = 8,height = 8)
par(mar=c(2,2,1,1))
corrplot(as.matrix(res$r), p.mat=as.matrix(res$P), sig.level = 0.05, insig="label_sig", diag = FALSE, type = 'lower')#,upper="ellipse",title="Corrplot")
dev.off()

##  regress seqPC1-4

ToKeep = c("seqPC1", "seqPC2", "seqPC3", "seqPC4")
designvars <- AlignBioCov[,ToKeep,drop=F]
expression <- "lm(thisExpr ~ seqPC1 + seqPC2 + seqPC3 + seqPC4, data = designmatrix)"
Method <- "lm"

LIST <- lm.fun(Ex,expression,designvars)
regEx <- LIST[[1]]


res <- cal_cov_ex_corr(regEx, AlignBioCov)
pdf(file=paste0(outdir,'/cov_ex_seqPC1-4_corrplot.pdf'),width = 8,height = 8)
par(mar=c(2,2,1,1))
corrplot(as.matrix(res$r), p.mat=as.matrix(res$P), sig.level = 0.05, insig="label_sig", diag = FALSE, type = 'lower')#,upper="ellipse",title="Corrplot")
dev.off()

Cov <- cbind(Cov, AlignBioCov%>%select(seqPC1, seqPC2, seqPC3, seqPC4))
save(Ex, geneAnno, Cov, file=file.path(outdir,'/preprocess.rda'))

