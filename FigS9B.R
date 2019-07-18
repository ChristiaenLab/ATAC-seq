library(edgeR)
library(DBI)
source('data/dirfns.R')
source('data/sqlfns.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
gene.names <- getGeneNames(con)
rna <- getRnaDat(con)

Foxf_LacZ <- rna$FoxF10hpf_LacZ10hpf

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
