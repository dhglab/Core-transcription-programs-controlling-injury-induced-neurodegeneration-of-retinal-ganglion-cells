rm(list=ls())

source("./code/00_utils.R")

### 1. DEGs of each TF ####
outdir <- './output'
load(file.path(outdir, '/DEX.rda'))

DEX.p <- DEX[,grep('FDR',colnames(DEX))]
DEX.logFC <- DEX[,grep('logFC',colnames(DEX))]

geneList <- geneAnno$ensembl_gene_id
dex.gene.list <- cl.list <- NULL 
for(tf in tfs){
  
  idx1 <- DEX.logFC[,paste0(tf,'.logFC')] > 0 & DEX.p[,paste0(tf,'.FDR')] < 0.1
  idx2 <- DEX.logFC[,paste0(tf,'.logFC')] < 0 & DEX.p[,paste0(tf,'.FDR')] < 0.1
  
  dex1 <- geneList[idx1] %>% na.omit() %>% as.character();cl1 <- rep('Up',length(dex1));names(cl1) <- dex1
  dex2 <- geneList[idx2] %>% na.omit() %>% as.character();cl2 <- rep('Down',length(dex2));names(cl2) <- dex2
  
  dex.gene <- c(dex1,dex2)
  cl <- c(cl1,cl2)
  
  
  dex.gene.list[[tf]] <- dex.gene
  cl.list[[tf]] <- cl
  
  
}

### 2. footprint of each TF ####

ft.gene.list <- NULL
for(tf in c("ATF3", "ATF4", "CEBPG","DDIT3")){
  ft <- read.table(file = paste0('./data/',tf,'.merged.filtered.geneBody.txt'), sep = '\t',header = T,stringsAsFactors=FALSE,quote = "")
  idx <- !is.na(ft$SYMBOL) #remove NA symbols
  ft <- ft[idx,]
  
  ft.gene.list[[tf]] <- unique(ft$"geneId")
  
}
names(ft.gene.list)[4] <- 'CHOP'


### 3.  overlap of each tf's dex and ft ####
target.list <- logfc.list <-  NULL 

for(tf in tfs){
  
  this.target  <- intersect(ft.gene.list[[tf]], dex.gene.list[[tf]]) # intersect analysis 
  this.target.gene <- geneAnno$gene[match(this.target, geneAnno$ensembl_gene_id)]
  
  write.table(DEX[this.target.gene,paste0(tf,'.logFC'), drop=F], 
              file = paste0(outdir,'/', tf, '_ft_deg.txt'), sep = '\t', quote = F,col.names = F)
  
  this.logFC <- DEX[this.target.gene,paste0(tf,'.logFC')]
  names(this.logFC) <- this.target.gene
  
  logfc.list[[tf]] <- this.logFC
  target.list[[tf]] <- this.target
}


save(target.list, logfc.list,ft.gene.list,dex.gene.list, DEX, Cov, geneAnno, file=paste0(outdir, '/ft.dex.rda'))



### 4. genes unique or shared by ATF3/CHOP, ATF4/CEBPG and GO ####

for(select in c('up', 'down')){
  if(select == 'down'){
    this.logfc.list <- lapply(logfc.list, function(x){ x[x<0]})}else{
      this.logfc.list <- lapply(logfc.list, function(x){ x[x>0]})
    }
  
  
  ##union of atf3 and chop targets
  target1 <- c(names(this.logfc.list[['ATF3']]),names(this.logfc.list[['CHOP']])) %>% 
    unlist() %>% unique()
  ## union atf4 and cebpg targets
  target2 <- c(names(this.logfc.list[['ATF4']]),names(this.logfc.list[['CEBPG']])) %>% 
    unlist() %>% unique()
  
  ## shared by two sets
  target3 <-intersect(target1, target2)
  
  ## unique to ATF3 and CHOP
  target11 <- target1[is.na(match(target1,target3))]
  
  ## unique to ATF4 and CEBPG
  target22 <- target2[is.na(match(target2,target3))]
  
  this.target.list <- list('ATF3:CHOP' = target11, 'ATF4:CEBPG' = target22, 'Common'=target3)
  target <- plyr::ldply(this.target.list, cbind)
  colnames(target) <- c('term_id','genes')
  target <- target %>% mutate(direction = select)
  
  #for IPA analysis
  write.table(as.data.frame(target11), file = paste0(outdir, '/ATF3CHOP_unique_', select,'.txt'), sep = '\t', quote = F, row.names = F,col.names = F)
  write.table(as.data.frame(target22), file = paste0(outdir, '/ATF4CEBPG_unique_',select,'.txt'), sep = '\t', quote = F, row.names = F,col.names = F)
  write.table(as.data.frame(target3), file = paste0(outdir, '/commonby2_', select,'.txt'), sep = '\t', quote = F, row.names = F,col.names = F)
  write.table(as.data.frame(target), file = paste0(outdir, '/All_', select,'.txt'), sep = '\t', quote = F)
  
  

}

