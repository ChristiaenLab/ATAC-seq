# Fig. 2B, S9D
source('data/sqlfns.R')
source('data/dirfns.R')
source('data/DESeqFns.R')

library(DBI)

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

bulkGS <- getBulkRNA(con)
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
peakGeneAnnotation <- getAnnotation(con)
gene.names <- peakGeneAnnotation$gene.names
# Fig. 2B

atac <- getAtac(con)

mespDnSig <- sig.sub(atac$condition_mesp_dnFGFR_vs_control,.5)
timeSig <- sig.sub(atac$time_6hpf_vs_10hpf,.5)
foxfSig <- sig.sub(atac$condition_FoxF_KO_vs_control,.45)

scrnaGenePeak <- lapply(
  prime.denovo,
  mergeGenePeak,
  intersect(row.names(timeSig),row.names(mespDnSig)),
  con=con
)

dir.eps('mesp_dn_6_10')
plot(
  atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,
  atac$time_6hpf_vs_10hpf$log2FoldChange,
  pch=19,col='gray',xlim=c(-2.5,2.5),ylim=c(-4,4),
  cex=1.5,cex.axis=1.5
)

points(
  mespDnSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  pch=19,col='gray28',cex=1.5
)
points(
  mespDnSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  pch=19,col='forestgreen',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  pch=19,col='blue',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  pch=19,col='red',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  pch=1,col='blue',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  pch=1,col='red',cex=1.5
)
lines(c(0,0),c(-4,4))
lines(c(-2.5,2.5),c(0,0))
dev.off()
# Fig. 2B Spearman correlation
cor(
  atac$condition_mesp_dnFGFR_vs_control[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
    ],
  atac$time_6hpf_vs_10hpf[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
  ],
  method = 'spearman'
)

# Fig. S9D
scrnaGenePeak <- lapply(
  prime.denovo,
  mergeGenePeak,
  intersect(row.names(foxfSig),row.names(mespDnSig)),
  con=con
)
dir.eps('mesp_dn_foxf_ko')
plot(
  atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,
  atac$condition_FoxF_KO_vs_control$log2FoldChange,
  pch=19,col='gray',xlim=c(-2.5,2.5),ylim=c(-2,2),
  cex=1.5,cex.axis=1.5
)
points(
  mespDnSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  pch=19,col='gray28',cex=1.5
)
points(
  mespDnSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  pch=19,col='forestgreen',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  pch=19,col='blue',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  pch=19,col='red',cex=1.5
)
points(
  mespDnSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  pch=1,col='blue',cex=1.5
)

points(
  mespDnSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  pch=1,col='red',cex=1.5
)
lines(c(0,0),c(-2,2))
lines(c(-2.5,2.5),c(0,0))
dev.off()
# 
dir.eps('foxf_6_10')
plot(
  atac$condition_FoxF_KO_vs_control$log2FoldChange,
  atac$time_6hpf_vs_10hpf$log2FoldChange,
  pch=19,col='gray',xlim=c(-2.5,2.5),ylim=c(-4,4),
  cex=1.5,cex.axis=1.5
)

points(
  foxfSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  pch=19,col='gray28',cex=1.5
)
points(
  foxfSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  pch=19,col='forestgreen',cex=1.5
)

points(
  foxfSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  pch=19,col='blue',cex=1.5
)

points(
  foxfSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  pch=19,col='red',cex=1.5
)

points(
  foxfSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  pch=1,col='blue',cex=1.5
)

points(
  foxfSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  pch=1,col='red',cex=1.5
)
lines(c(0,0),c(-4,4))
lines(c(-2.5,2.5),c(0,0))
dev.off()
# Fig. S9D Spearman correlation
cor(
  atac$condition_mesp_dnFGFR_vs_control[
    peaksets$mespDep,'log2FoldChange'
    ],
  atac$condition_FoxF_KO_vs_control[
    peaksets$mespDep,'log2FoldChange'
  ],
  method = 'spearman'
)

cor(
  atac$condition_FoxF_KO_vs_control[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
    ],
  atac$time_6hpf_vs_10hpf[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
  ],
  method = 'spearman'
)
