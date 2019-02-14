library(edgeR)
library(DBI)
source('data/dirfns.R')
source('data/sqlfns.R')
source('~/Dropbox/ciona/code/plotfns.R')

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

Foxf_LacZ <- res$FoxF10hpf_LacZ10hpf

sel <- abs(Foxf_LacZ$log2FoldChange)>.75&Foxf_LacZ$padj<.05
heart <- intersect(row.names(Foxf_LacZ)[sel],scrna$Cardiac)
asm <- intersect(row.names(Foxf_LacZ)[sel],scrna$ASM)
tvcp <- intersect(row.names(Foxf_LacZ)[sel],scrna$TVCP)

# Fig. S9B
gene.names[heart,]
gene.names[asm,]
gene.names[tvcp,]

dir.eps('FoxF_LacZ_volcano')
plot(
  Foxf_LacZ$log2FoldChange[Foxf_LacZ$p.value<.10],
  -log10(Foxf_LacZ$padj[Foxf_LacZ$p.value<.10]),
  col='gray',pch=19,xlim=c(-4,4),cex=1.5
)
points(Foxf_LacZ[heart,"log2FoldChange"],-log10(Foxf_LacZ[heart,"padj"]),pch=19,col='red',cex=1.5)
points(Foxf_LacZ[asm,"log2FoldChange"],-log10(Foxf_LacZ[asm,"padj"]),pch=19,col='blue',cex=1.5)
points(Foxf_LacZ[tvcp,"log2FoldChange"],-log10(Foxf_LacZ[tvcp,"padj"]),pch=19,col='forestgreen',cex=1.5)
dev.off()


dbWriteGenes(con,'foxf_rnaseq',melt.rename(res,'comparison'))