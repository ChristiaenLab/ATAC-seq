library(DESeq2)
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
expDesign$condtime <- sub('.*lacz','lacz',expDesign$condtime)

dds = DESeqDataSetFromMatrix(
  countData=dat, 
  colData = expDesign, 
  design = ~condtime
)
dds = DESeq(dds)

comp <- list(
  c('handrdnfgfr12hpf','lacz12hpf'),c('foxfcamras12hpf','lacz12hpf'),
  c('handrdnfgfr15hpf','lacz15hpf'),c('foxfcamras15hpf','lacz15hpf'),
  c('handrdnfgfr18hpf','lacz18hpf'),c('foxfcamras18hpf','lacz18hpf'),
  c('handrdnfgfr20hpf','lacz20hpf'),c('foxfcamras20hpf','lacz20hpf'),
  c('lacz12hpf','lacz15hpf'),c('lacz12hpf','lacz20hpf'),
  c('handrdnfgfr18hpf','foxfcamras18hpf')
)
res <- lapply(comp,function(x) results(
  dds,c('condtime',x[1],x[2]),format = "DataFrame",alpha = .05
))
names(res) <- sapply(comp,function(x) paste(c("condtime",x),collapse = '_'))

dbWriteGenes(con,'handr_rnaseq',melt.rename(res,'comparison'))

up12 <- row.names(subset(res[[9]],padj<.05&log2FoldChange>1))
down12 <- row.names(subset(res[[9]],padj<.05&log2FoldChange< -1))
downFgfr15 <- row.names(subset(res[[3]],padj<.05&log2FoldChange< -1))
upFgfr15 <- row.names(subset(res[[3]],padj<.05&log2FoldChange>1))
downMek15 <- row.names(subset(res[[4]],padj<.05&log2FoldChange< -1))
upMek15 <- row.names(subset(res[[4]],padj<.05&log2FoldChange>1))
downFgfr12 <- row.names(subset(res[[1]],padj<.05&log2FoldChange< -1))
upFgfr12 <- row.names(subset(res[[1]],padj<.05&log2FoldChange>1))
downMek12 <- row.names(subset(res[[2]],padj<.05&log2FoldChange< -1))
upMek12 <- row.names(subset(res[[2]],padj<.05&log2FoldChange>1))

writeLines(setdiff(intersect(down12,downFgfr15),downFgfr12),'up15downFgfr.txt')
writeLines(setdiff(intersect(down12,downMek15),downMek12),'up15downMek.txt')
