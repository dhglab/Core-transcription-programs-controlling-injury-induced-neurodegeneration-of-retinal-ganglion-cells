
# differential analysis 

source("./code/00_utils.R")

load('./output/preprocess.rda')
outdir <- './output'

# remove outliers identified from 1b_qc_outlier.R
toDel <- c('ATF4.3d.6','ATF4.3d.7','ATF4.3d.5', 'CEBPG.3d.7', "non.3d.9","non.3d.11")
idx <- is.na(match(rownames(Cov), toDel))
Cov <- Cov[idx,]
Ex <- Ex[, rownames(Cov)]

save(Ex, Cov, geneAnno, file=paste0(outdir, '/filterEx.rda'))

##edgeR to filter lowly-expressed genes and normalize sequencing depth

y <- edgeR.process(Cov, Ex, geneAnno) 

normEx <- cpm(y)
geneAnno <- y$gene
Cov <- y$samples
rownames(normEx) <- geneAnno$gene

save(normEx, geneAnno, Cov,  file=paste0(outdir, '/normEx.rda'))

# differential gene analysis


design <- ~0+seqPC1+seqPC2+seqPC3+seqPC4+Target
mod <- model.matrix(as.formula(design), data=y$samples)
colnames(mod) <- gsub('Target','',colnames(mod))

y <- estimateDisp(y, mod, robust=TRUE)

pdf(paste0(outdir, '/dispersion.pdf'), width = 5,height = 5)
plotBCV(y)
dev.off()

fit <- glmQLFit(y, mod)


DEsets <- fitsets <- list()

for(this.tf in tfs){
  this.con <- paste0(this.tf, '.3d_vs_non.3d=', this.tf, '-non')
  fitsets[[this.tf]] <- glmQLFTest(fit, contrast=makeContrasts(this.con, levels=mod))
  
}

DEsets <- lapply(fitsets, function(i){topTags(i,n = Inf, adjust.method = 'fdr',sort.by = 'none')})
names(DEsets) <- names(fitsets)<- tfs

DEGlist <- lapply(DEsets, function(x){x@.Data[[1]][,-c(1:5)]})
DEX <- do.call(cbind, DEGlist)
info(DEX)
rownames(DEX) <- geneAnno$gene

geneName <- geneAnno$gene
genes <- bitr(geneName,fromType = "SYMBOL",toType = c("ENSEMBL"),OrgDb = org.Mm.eg.db)
geneAnno$ensembl_gene_id <- genes$ENSEMBL[match(geneAnno$gene, genes$SYMBOL)]

save(DEsets, fitsets, DEX,geneAnno, Cov,  file = paste0(outdir, '/DEX.rda'))




