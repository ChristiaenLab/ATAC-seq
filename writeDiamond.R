source('data/dbDiamond.R')
source('data/sqlfns.R')
source('data/dirfns.R')
source('scatterplotfns.R')
library(DBI)
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
bulkGS <- getBulkRNA(con)
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
peakGeneAnnotation <- getAnnotation(con)
gene.names <- mcols(peakGeneAnnotation$genes,use.names = T)

rna <- getRnaDat(con)[c(3:5,2)]
rna$gfp.lacz <- rna$FoxF10hpf_LacZ10hpf
rna$gfp.lacz$log2FoldChange <- 0
rna$gfp.lacz$padj <- 1
atacdat <- getAtacLib(con,c(
  "condition_handr_dnFGFR_vs_control","condition_handr_MekMut_vs_control",
  "condition_mesp_dnFGFR_vs_control","condition_FoxF_KO_vs_control",
  "tissue_B7.5_vs_mesenchyme"
))

prime.denovo <- scrna[c(
  'primedCardiac','primedASM','denovoCardiac','denovoASM','TVCP','STVC','ATM','mesenchyme'
)]
ebf <- scrna[c('ebfActivated','ebfInhibited')]
cols <- c(
  primedCardiac='red',primedASM='blue',
  denovoCardiac='hotpink',denovoASM='skyblue',
  TVCP='forestgreen',STVC='orange',ATM='gray28',
  ebfActivated='lightgreen',ebfInhibited='purple'
)

tmp <- setNames(reshape2::melt(
  mapply(
    function(genes,peaks) setdiff(genes,mergeGenePeak(con,genes,peaks)$GeneID),
    genes=bulkGS[c("MAPK18activated","MAPK18inhibited")],
    peaks=peaksets[c('asmAcc','heartAcc')]
  )
),c("GeneID","gene_set"))
tmp$UniqueNAME <- gene.names[as.character(tmp$GeneID),"UniqueNAME"]
dir.tab(tmp,'MAPK18noDA',row.names=F)

dbDiamondplot(
  con,rna,atac,
  prime.denovo[c('denovoCardiac','denovoASM')],
  list(
    peaksets$heartAcc,
    peaksets$asmAcc,
    intersect(peaksets$heartAcc,peaksets$open18),
    intersect(peaksets$asmAcc,peaksets$open18)
  ),
  c('red','blue','hotpink','skyblue'),
  list(prime.denovo),cols,'mapk18prime.denovo',T
)

dbDiamondplot(
  con,rna,atac,
  sapply(prime.denovo[c('denovoCardiac','denovoASM')],intersect,union(
    row.names(sig.sub(rna$condtime_foxfcamras18hpf_handrlacz18hpf)),
    row.names(sig.sub(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf))
  )),
  list(heartAcc,asmAcc),c('red','blue'),list(prime.denovo),cols,'mapk18denovo',F
)

dbDiamondplot(
  con,rna,atac,
  bulkGS[c(
    "MAPK18activated","MAPK18inhibited"
  )],
  # list(
  #   MAPK18activated=intersect(bulkGS$MAPK18activated,scrna$denovoASM),
  #   MAPK18inhibited=intersect(bulkGS$MAPK18inhibited,scrna$denovoCardiac)
  # ),
  list(heartAcc,asmAcc),c('red','blue'),list(prime.denovo),cols,'mapk18',F
)
dbDiamondplot(
  con,rna,atac,
  sapply(
    bulkGS[c(
      "FoxFactivated","FoxFinhibited"
    )],
    function(x) mergeGenePeak(
      con,x,peaksets$closedFoxf
    )$GeneID
  ),
  list(peaksets$tvcAcc,intersect(peaksets$tvcAcc,peaksets$closed6)),
  c('forestgreen','brown'),
  list(prime.denovo),cols,'foxf',T
)

quant.lfc <- c(
  quantile(atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,.0075),
  quantile(atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,.9925)
)

mapk10quant <- lapply(
  bulkGS[c(
    "MAPK10activated","MAPK10inhibited"
  )],
  function(x) unique(mergeGenePeak(
    con,
    x,
    # intersect(x,bulkGS$downreg6hpf),
    row.names(sig.sub(atac$condition_mesp_dnFGFR_vs_control,quant.lfc))
  )$GeneID)
)

dbDiamondplot(
  con,rna,atac,
  mapk10quant,
  list(tvcAcc,atmAcc,heartAcc,asmAcc),
  c('forestgreen','gray28','red','blue'),
  list(prime.denovo),cols,'mapk10'
)
dbDiamondplot(
  con,rna,atac,
  lapply(mapk10quant,intersect,bulkGS$downreg6hpf),
  list(tvcAcc,atmAcc,heartAcc,asmAcc),
  c('forestgreen','gray28','red','blue'),
  list(prime.denovo),cols,'mapk10downreg6'
)

dbDiamondplot(
  con,rna,atac,
  list(
    MAPKinhibited=union(bulkGS$MAPK18inhibited[
      order(rna$condtime_foxfcamras18hpf_handrlacz18hpf[
        bulkGS$MAPK18inhibited,'log2FoldChange'
      ])[1:50]
    ],intersect(prime.denovo$TVCP,bulkGS$MAPK18inhibited)),
    MAPKactivated=union(bulkGS$MAPK18activated[
      order(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf[
        bulkGS$MAPK18activated,'log2FoldChange'
      ])[1:50]
    ],intersect(prime.denovo$STVC,bulkGS$MAPK18activated))
  ),
  list(heartAcc,asmAcc),
  c('red','blue'),
  list(prime.denovo),cols,'mapk18top50',F
)

tmp <- mergeGenePeak(con,bulkGS$MAPK10activated,peaksets$tvcAcc)
tmp$UniqueNAME <- gene.names[tmp$GeneID,"UniqueNAME"]
tmp <- merge(tmp,dbReadTable(con,'peakfeature')[,1:4])
dir.tab(tmp,'mapk10activated_tvcAcc',row.names=F)

scrnaPeaks <- sapply(prime.denovo,mergeGenePeak,con=con)["PeakID",]

dir.eps('mesp_dn_6_10')
plot(
  atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,
  atac$time_6hpf_vs_10hpf$log2FoldChange,
  pch=19,col='gray',xlim=c(-2.5,2.5),ylim=c(-4,4)
)
# mapply(
#   function(x,col,pch) points(
#     atac$condition_mesp_dnFGFR_vs_control[x,"log2FoldChange"],
#     atac$time_6hpf_vs_10hpf[x,"log2FoldChange"],
#     col=col,pch=pch
#   ),
#   x=scrnaPeaks[c(
#     "ATM","TVCP","primedCardiac","primedASM","denovoCardiac","denovoASM"
#   )],
#   col=c('gray28','forestgreen','red','blue','red','blue'),
#   pch=c(19,19,19,19,1,1)
# )
mespDnSig <- sig.sub(atac$condition_mesp_dnFGFR_vs_control,.5)
timeSig <- sig.sub(atac$time_6hpf_vs_10hpf,.5)
foxfSig <- sig.sub(atac$condition_FoxF_KO_vs_control,.45)
points(
  mespDnSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  pch=19,col='gray28'
)
points(
  mespDnSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  pch=19,col='forestgreen'
)

points(
  mespDnSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  pch=19,col='blue'
)

points(
  mespDnSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  pch=19,col='red'
)

points(
  mespDnSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  pch=1,col='blue'
)

points(
  mespDnSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  timeSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  pch=1,col='red'
)
lines(c(0,0),c(-4,4))
lines(c(-2.5,2.5),c(0,0))
# points(
#   atac$condition_mesp_dnFGFR_vs_control[intersect(peaksets$timeDep,peaksets$mespDep),"log2FoldChange"],
#   atac$time_6hpf_vs_10hpf[intersect(peaksets$timeDep,peaksets$mespDep),"log2FoldChange"],
#   col='brown',pch=19
# )
dev.off()
cor(
  atac$condition_mesp_dnFGFR_vs_control[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
    ],
  atac$time_6hpf_vs_10hpf[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
  ],
  method = 'spearman'
)

dir.eps('mesp_dn_foxf_ko')
plot(
  atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,
  atac$condition_FoxF_KO_vs_control$log2FoldChange,
  pch=19,col='gray',xlim=c(-2.5,2.5),ylim=c(-2,2)
)
# mapply(
#   function(x,col,pch) points(
#     atac$condition_mesp_dnFGFR_vs_control[x,"log2FoldChange"],
#     atac$time_6hpf_vs_10hpf[x,"log2FoldChange"],
#     col=col,pch=pch
#   ),
#   x=scrnaPeaks[c(
#     "ATM","TVCP","primedCardiac","primedASM","denovoCardiac","denovoASM"
#   )],
#   col=c('gray28','forestgreen','red','blue','red','blue'),
#   pch=c(19,19,19,19,1,1)
# )
points(
  mespDnSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$ATM$PeakID,"log2FoldChange"],
  pch=19,col='gray28'
)
points(
  mespDnSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$TVCP$PeakID,"log2FoldChange"],
  pch=19,col='forestgreen'
)

points(
  mespDnSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$primedASM$PeakID,"log2FoldChange"],
  pch=19,col='blue'
)

points(
  mespDnSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$primedCardiac$PeakID,"log2FoldChange"],
  pch=19,col='red'
)
points(
  mespDnSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$denovoASM$PeakID,"log2FoldChange"],
  pch=1,col='blue'
)

points(
  mespDnSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  foxfSig[scrnaGenePeak$denovoCardiac$PeakID,"log2FoldChange"],
  pch=1,col='red'
)
lines(c(0,0),c(-2,2))
lines(c(-2.5,2.5),c(0,0))
dev.off()

cor(
  atac$condition_mesp_dnFGFR_vs_control[
    peaksets$mespDep,'log2FoldChange'
    ],
  atac$condition_FoxF_KO_vs_control[
    peaksets$mespDep,'log2FoldChange'
  ],
  method = 'spearman'
)

tmp <- mergeGenePeak(
  con,
  row.names(sig.sub(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf)),
  row.names(sig.sub(atac$condition_handr_dnFGFR_vs_control,.5))
)
cor(
  rna$condtime_handrdnfgfr18hpf_handrlacz18hpf[tmp$GeneID,'log2FoldChange'],
  atac$condition_handr_dnFGFR_vs_control[tmp$PeakID,'log2FoldChange'],
  method='spearman'
)

tmp <- mergeGenePeak(
  con,
  row.names(sig.sub(rna$FoxF10hpf_LacZ10hpf,.75)),
  row.names(sig.sub(atac$condition_FoxF_KO_vs_control,.45))
)
cor(
  rna$FoxF10hpf_LacZ10hpf[tmp$GeneID,'log2FoldChange'],
  atac$condition_FoxF_KO_vs_control[tmp$PeakID,'log2FoldChange'],
  method='spearman'
)
