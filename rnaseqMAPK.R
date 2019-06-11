library(DESeq2)
library(DBI)
source("data/dirfns.R")
source("data/sqlfns.R")
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

counts <- lrtab('reads_count',read.delim,row.names=1)
dat <- as.data.frame(
  sapply(counts,'[',,'accepted_hits.bam',drop=F),
  row.names = row.names(counts[[1]])
)
names(dat) <- sub('\\..*','',names(dat))

expDesign <- data.frame(
  condtime=sub('[0-9]+$','',names(dat)),
  row.names = names(dat)
)
# expDesign$experiment <- ifelse(grepl('^handr',expDesign$condtime),'handr','foxf')
# expDesign$condtime <- sub('.*lacz','lacz',expDesign$condtime)

dds <- DESeqDataSetFromMatrix(
  countData=dat, 
  colData = expDesign, 
  design = ~condtime
)
dds <- DESeq(dds)

comp <- list(
  c('handrdnfgfr12hpf','handrlacz12hpf'),c('foxfcamras12hpf','handrlacz12hpf'),
  c('handrdnfgfr15hpf','handrlacz15hpf'),c('foxfcamras15hpf','handrlacz15hpf'),
  c('handrdnfgfr18hpf','handrlacz18hpf'),c('foxfcamras18hpf','handrlacz18hpf'),
  c('handrdnfgfr20hpf','handrlacz20hpf'),c('foxfcamras20hpf','handrlacz20hpf'),
  c('handrlacz12hpf','handrlacz15hpf'),c('handrlacz12hpf','handrlacz20hpf'),
  c('foxfcamras18hpf','handrdnfgfr18hpf')
)
res <- lapply(comp,function(x) results(
  dds,c('condtime',x[1],x[2]),format = "DataFrame",alpha = .05
))
res <- lapply(res,function(x) cbind(GeneID=row.names(x),x))
names(res) <- sapply(comp,function(x) paste(c("condtime",x),collapse = '_'))

dbWriteGenes(con,'handr_rnaseq',melt.rename(res,'comparison'))

up12 <- subset(res[[9]],padj<.05&log2FoldChange>1)$GeneID
down12 <- subset(res[[9]],padj<.05&log2FoldChange< -1)$GeneID
downFgfr15 <- subset(res[[3]],padj<.05&log2FoldChange< -1)$GeneID
upFgfr15 <- subset(res[[3]],padj<.05&log2FoldChange>1)$GeneID
downMek15 <- subset(res[[4]],padj<.05&log2FoldChange< -1)$GeneID
upMek15 <- subset(res[[4]],padj<.05&log2FoldChange>1)$GeneID
downFgfr12 <- subset(res[[1]],padj<.05&log2FoldChange< -1)$GeneID
upFgfr12 <- subset(res[[1]],padj<.05&log2FoldChange>1)$GeneID
downMek12 <- subset(res[[2]],padj<.05&log2FoldChange< -1)$GeneID
upMek12 <- subset(res[[2]],padj<.05&log2FoldChange>1)$GeneID

dir.tab(mergeGenePeak2(con,setdiff(intersect(down12,downFgfr15),downFgfr12),peaksets$tvcAcc),'up15downFgfrTVCacc',row.names=F)

dir.tab(merge(
  mergeGenePeak2(con,setdiff(intersect(down12,upFgfr15),upFgfr12),peaksets$tvcAcc),
  ann$features,
  by.x="PeakID",
  by.y="row.names"
),'up15upFgfrTVCacc',row.names=F)

dir.tab(merge(
  mergeGenePeak2(con,setdiff(intersect(down12,downMek15),downMek12),peaksets$tvcAcc),
  ann$features,
  by.x="PeakID",
  by.y="row.names"
),'up15downMekTVCacc',row.names=F)
