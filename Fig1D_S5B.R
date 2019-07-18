library(fgsea)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(lattice)

source('data/dirfns.R')

gsea <- Reduce(
  function(y,atac){
    load(paste(atac,'fgseaRes.Rdata',sep='/'))
    row.names(fgseaRes$fgsea) <- fgseaRes$fgsea$pathway
    y[[atac]] <- fgseaRes$fgsea
    return(y)
  },
  list.files('2018-12-14/gsea/peaksets/',full.names = T),
  init = NULL
)
names(gsea) <- sub('\\.NES','',sub('.*/','',names(gsea)))
nes <- do.call(cbind,sapply(gsea,'[',,'NES'))
gsea.fdr <- do.call(cbind,sapply(gsea,'[',,'padj'))
row.names(nes) <- gsea[[1]]$pathway
row.names(gsea.fdr) <- row.names(nes)

bar <- barchart(
  value~Var1|Var2,
  melt(nes[
    c("ATM","TVCP","primedCardiac","primedASM","denovoCardiac","denovoASM"),
  ]),groups=melt(gsea.fdr[
    c("ATM","TVCP","primedCardiac","primedASM","denovoCardiac","denovoASM"),
  ]<.05)$value,
  origin=0,box.ratio = 8,
  col=c('gray','gray28'),
  scales=list(x=list(rot=90))
)
dir.eps('Fig1D',height=24,width=12)
plot(bar)
dev.off()

nes[gsea.fdr>.05] <- NaN

dir.eps('FigS5B')
Heatmap(
  nes[c(
    "MAPK18inhibited","MAPK18activated","MAPK10inhibited","MAPK10activated",
    "FHP14","FHP","SHP","STVC","TVCP","ATM",
    "primedASM","denovoASM","primedCardiac","denovoCardiac"
  ),c(
    "condition_mesp_dnFGFR_vs_control.NES",
    "time_6hpf_vs_10hpf.NES","time_10hpf_vs_15hpf.NES","time_10hpf_vs_18hpf.NES",
    "tissue_B7.5_vs_mesenchyme.NES"
  )],
  cluster_columns = F,
  cluster_rows = F
)
dev.off()