library(edgeR)
library(DBI)
source('data/dirfns.R')
source('data/sqlfns.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
gene.names <- getGeneNames(con)

counts <- lrtab('star',pattern = 'count.txt',row.names=1)
dat <- do.call(cbind,counts)
names(dat) <- sub('_count','',names(counts))
expDesign <- data.frame(
  condition=c(rep('FoxF10hpf',2),rep("LacZ10hpf",2),rep("Ngn10hpf",2)),
  date=sub('.*rnaseq_','',names(dat)),
  row.names = names(dat)
)
filt <- dat[-grep('_',row.names(dat)),]

comp <- matrix(c('FoxF10hpf',"LacZ10hpf","FoxF10hpf","Ngn10hpf"),2)

design <- model.matrix(~ 0+condition+date,expDesign)
dge <- DGEList(filt,group=expDesign$condition,samples = expDesign)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge,design)

res <- apply(comp,2,function(x) setNames(
  exactTest(dge,c(x[2],x[1]))$table,
  c('log2FoldChange','logCPM','p.value')
))
res <- lapply(res,function(x) cbind(GeneID=row.names(x),x,padj=p.adjust(x$p.value)))

names(res) <- apply(comp,2,paste0,collapse='_')

dbWriteGenes(con,'foxf_rnaseq',melt.rename(res,'comparison'))

