##++++++++ TOBIAS footprint process Step1 +++++++++##
## 1. filter footprint
##    -logFC > 0 for activators and logFC <0 for repressors
##    - at least bound in one condition
## 2. Merge TFBS from multiple sources, hocomocco, jaspar
##   
## 3. Extend TFBS from summit +-100bp to extract counts for comparison
##++++++++++++++++++++++++++++++++++++++++++++++++++##

#module load R/3.6.1
library(GenomicFeatures)
library(ChIPseeker)


dir='/u/home/y/ycheng41/nobackup-dhg/He_RGC/ATAC/CTCF/nonTarget/TOBIAS/BINDetect'
#suppressor.tfs <- c('CTCF','REST')
#activator.tfs <- c('ATF3','ATF4','CEBPG','DDIT3','CEBPA')
core.tfs <- c('ATF3','ATF4','CEBPG','DDIT3')

##1. filter footprints
#   -logFC > 0 for activators and logFC <0 for repressors
#   - at least bound in one condition

#BINDetec output files are TFBS in consensus peaks
# select peaks that has footprint (bound) in either day 1 or 3, but no bound in 0
indir=paste0(dir,'/overview')
outdir=paste0(dir,'/annotated_tfbs') #outdir=paste0(dir,'/annotated_tfbs/filter') #for 

if(!dir.exists(outdir)){dir.create(outdir)}

setwd(indir)
samplefiles <- list.files(indir, pattern="_overview.txt", full.names=F) #$means the end of filename
idx <- grep(paste0(core.tfs,collapse = '|'),samplefiles, ignore.case = T )
samplefiles <- samplefiles[idx]

name.prefix <- basename(samplefiles)
name.prefix <- gsub('var.2','',name.prefix)
name.prefix <- gsub("\\..*", "", name.prefix)
name.prefix <- toupper(name.prefix)

name.tf <- gsub("\\_.*", "", name.prefix)



for(i in 1:length(samplefiles)){
  ft <- read.table(samplefiles[i], header = T, sep = '\t')
  colnames(ft) <- gsub(".unique.sorted.rmDup.merge.footprint","",colnames(ft))
  
  colnames(ft) <- gsub("_b_","_",colnames(ft))
  
  ft$nonTarget_day0_nonTarget_day3_log2fc <- -ft$nonTarget_day0_nonTarget_day3_log2fc
  colnames(ft)[18] <- 'nonTarget_day3_nonTarget_day0_log2fc'
  
  #any of injured has more activty than control
  idx <- ft$nonTarget_day1_nonTarget_day0_log2fc >0 & ft$nonTarget_day3_nonTarget_day0_log2fc >0
  ft <- ft[idx,]
  
  #unbound in control but bound in injury
  #idx2 <- ft$nonTarget_day0_bound == 0 & rowSums(ft[,grep('bound',colnames(ft))]) >0
  idx2 <- rowSums(ft[,grep('bound',colnames(ft))]) >0
  ft <- ft[idx2,]
  
  write.table(ft, file = paste0(outdir,'/',name.prefix[i],'_filtered.bed'),
              col.names = F, row.names = F, sep = '\t', quote = F)

}

##2. merge footprints from hocomocco and jaspar and average the site scores 

setwd(outdir)
system (
  module add bedtools
  #for tf in CTCF REST SOX11 STAT3 JUN KLF6 FOSL2 BACH1 BHLHE41 ETS2 POU2F1 JUND ETV6 STAT2 BCL6 JUND JUNB KLF9 RELA IRF9 CREB1 HIF1A MYC; do
  for tf in ATF3 ATF4 CEBPG DDIT3; do
  mergedPeak=$tf'.merged.filtered.bed' 
  
  cat $tf* | sort -k1,1 -k2,2n | bedtools merge -i - -d -1 \
  -c 12,11,13,15,14,16,17,18 -o mean,mean,mean,max,max,max,mean,mean > $mergedPeak
  
  done
)


##3. Extend TFBS from summit +-100bp to extract featureCounts
#To avoid biases from TFBS clustering on the same ATAC-seq peak, choose the peak with the highest average counts across all samples
samplefiles <- list.files(outdir, pattern=".merged.filtered.bed", full.names=T) #$means the end of filename
idx <- grep(paste0(core.tfs,collapse = '|'),samplefiles, ignore.case = T )
samplefiles <- samplefiles[idx]
name.tf <- gsub("\\..*", "", basename(samplefiles))
names(samplefiles) <- name.tf

for(tf in name.tf){
  ft <- readPeakFile(samplefiles[tf])
  ft <- resize(ft, width = 200, fix = "center", ignore.strand=T) #extend 100bp each side from the summit
  
  ft <- as.data.frame(ft)
  #ft <- ft[,-c(4,5)] #remove strand and width information
  write.table(ft, file =paste0(outdir,'/',tf,'.merged.filtered.extend.bed'),
              col.names = F, row.names = F, sep = '\t', quote = F)
  
}

##4. get protmoer and distal regions of TFBS

#a. prepare TF-geneID annotation- use .mtf files in rgtdata folder
jaspar <- read.table(file="~/rgtdata/motifs/jaspar_vertebrates.mtf", sep = '\t',header = F) 
hocomoco <- read.table(file="~/rgtdata/motifs/hocomoco.mtf", sep = '\t',header = F)

TF <- unique(c(tolower(as.character(jaspar$V4)), tolower(as.character(hocomoco$V1))))
TF <- tools::toTitleCase(TF)
#convert uniport ID to ensembl 

TF_ID <- bitr(TF, fromType = "SYMBOL",toType = c("ENSEMBL","SYMBOL"),OrgDb = org.Mm.eg.db)
TF_ID$SYMBOL <- toupper(TF_ID$SYMBOL)
write.table(TF_ID, file = paste0(dir,'/motif2gene_mapping.txt'),
            row.names = F,col.names = F, sep = '\t',quote = F)

#b. get promoter regions for intput as create network - run in terminal

#prepare txdb from genecode.vM24 

gencode.vM24 <- "~/project-geschwind/gencode.vM24/gencode.vM24.annotation.gtf.gz"
txdb<- makeTxDbFromGFF(gencode.vM24, format="gtf", organism = "Mus musculus", 
                       chrominfo = seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene))


saveDb(txdb, '~/project-geschwind/gencode.vM24/txdb.gencode24.mm10')
#txdb = loadDb(file = '~/project-geschwind/gencode.vM24/txdb.gencode24.mm10')

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, annoDb='org.Mm.eg.db',
                       tssRegion=c(-2000, 2000), verbose=FALSE)

names(peakAnnoList) <- name.tf
name.tf <- factor(name.tf, levels = c('ATF3','ATF4','CEBPG','DDIT3'))
peakAnnoList <- peakAnnoList[order(name.tf)]
save(peakAnnoList,file = paste0(outdir,'/annotated.footprint.rda'))

pdf(file=paste0(outdir,"/footprintAnnoBar.pdf"), width = 4, height = 3.5)
plotAnnoBar(peakAnnoList)
dev.off()


pdf(file=paste0(outdir,"/footprintDis.pdf"), width = 4, height = 3.5)
plotDistToTSS(peakAnnoList,title="")
dev.off()

for(tf in names(peakAnnoList)){
  
  ft <- peakAnnoList[[tf]]
  ft <- as.data.frame(ft)
  
  colnames(ft)[6:8] <- paste0(c('nonTarget_day0','nonTarget_day1','nonTarget_day3'),'.footprint_score')
  colnames(ft)[9:11] <- paste0(c('nonTarget_day0','nonTarget_day1','nonTarget_day3'),'.bound')
  colnames(ft)[12:13] <- paste0(c('nonTarget_Day1vsDay0','nonTarget_Day3vsDay0'),'.logFC')
  ft$SYMBOL <- toupper(ft$SYMBOL)
  ft$geneId <- gsub("\\..*","", ft$geneId)
  ft$TFBS_name <- toupper(tf)
  
   write.table(ft, file = paste0(outdir,'/',tf,'.merged.filtered.txt'),
              sep = '\t', quote = F)
  
  #get promoter regions only
  idx <- grepl('Promoter', ft$annotation)
  ft_promoter <- ft[idx,]
  write.table(ft_promoter, file = paste0(outdir,'/',tf,'.merged.filtered.promoter.bed'),
              col.names = F, row.names = F, sep = '\t', quote = F)
  write.table(ft_promoter, file = paste0(outdir,'/',tf,'.merged.filtered.promoter.txt'),
              sep = '\t', quote = F)
  
  #get distal region
  idx <- grepl('Distal', ft$annotation)
  ft_distal <- ft[idx,]
  write.table(ft_distal, file = paste0(outdir,'/',tf,'.merged.filtered.distal.bed'),
              col.names = F, row.names = F, sep = '\t', quote = F)
  write.table(ft_distal, file = paste0(outdir,'/',tf,'.merged.filtered.distal.txt'),
              sep = '\t', quote = F)
  
  #get gene body
  
  ft_gb <- ft[!idx,]
  write.table(ft_gb, file = paste0(outdir,'/',tf,'.merged.filtered.geneBody.bed'),
              col.names = F, row.names = F, sep = '\t', quote = F)
  write.table(ft_gb, file = paste0(outdir,'/',tf,'.merged.filtered.geneBody.txt'),
              sep = '\t', quote = F)
  
  
}



