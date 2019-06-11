#Figs. 1EFG, 2CD,3AB, S5A, S10E, S17AB, S18B

source('data/dbDiamond.R')
source('data/sqlfns.R')
source('data/dirfns.R')
source('data/scatterplotfns.R')
source('data/DESeqFns.R')

library(DBI)

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

bulkGS <- getBulkRNA(con)
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
peakGeneAnnotation <- getAnnotation(con)
gene.names <- peakGeneAnnotation$gene.names

rna <- getRnaDat(con)[c(
  "condtime_handrdnfgfr18hpf_handrlacz18hpf","condtime_foxfcamras18hpf_handrlacz18hpf",
  "MA_dnFGFR_LacZ_10hpf","FoxF10hpf_LacZ10hpf"
)]
rna$gfp.lacz <- rna$FoxF10hpf_LacZ10hpf
rna$gfp.lacz$log2FoldChange <- NA
rna$gfp.lacz$padj <- 1
atac <- getAtacLib(con,c(
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

tmp <- mergeGenePeak(con,bulkGS$MAPK10activated,peaksets$tvcAcc)
tmp$UniqueNAME <- gene.names[tmp$GeneID,"UniqueNAME"]
tmp <- merge(tmp,dbReadTable(con,'peakfeature')[,1:4])
dir.tab(tmp,'mapk10activated_tvcAcc',row.names=F)

# Figs. 1G,5A
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
  peaksets[c("tvcAcc","atmAcc")],
  c('forestgreen','gray28'),
  list(prime.denovo),cols,'mapk10'
)

# Fig. 2D
dbDiamondplot(
  con,rna,atac,
  # sapply(
    bulkGS[c(
      "FoxFactivated","FoxFinhibited"
    )],
  #   function(x) mergeGenePeak(
  #     con,x,peaksets$closedFoxf
  #   )$GeneID
  # ),
  list(peaksets$tvcAcc,intersect(peaksets$tvcAcc,peaksets$closed6)),
  c('forestgreen','brown'),
  list(prime.denovo),cols,'foxf',
  gene.peak.intersect = F
)


# Figs. 3B, S18A 
dbDiamondplot(
  con,rna,atac,
  sapply(prime.denovo[c('denovoCardiac','denovoASM')],intersect,union(
    row.names(sig.sub(rna$condtime_foxfcamras18hpf_handrlacz18hpf)),
    row.names(sig.sub(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf))
  )),
  peaksets[c("heartAcc","asmAcc")],c('red','blue'),
  list(prime.denovo),cols,'mapk18denovo',
  gene.peak.intersect = F
)

# Fig. S17B
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
  peaksets[c("heartAcc","asmAcc")],
  c('red','blue'),
  list(prime.denovo),cols,'mapk18top50',
  gene.peak.intersect = F
)

# Figs. 1E, 2C, 3A, S17A
mapply(
  scatterRnaAtac,
  rna=rna[c(
    'condtime_foxfcamras18hpf_handrlacz18hpf',
    'condtime_handrdnfgfr18hpf_handrlacz18hpf',
    'FoxF10hpf_LacZ10hpf',
    'MA_dnFGFR_LacZ_10hpf'
  )],
  atac=atac[c(
    'condition_handr_MekMut_vs_control', 
    'condition_handr_dnFGFR_vs_control',
    'condition_FoxF_KO_vs_control',
    "condition_mesp_dnFGFR_vs_control"
  )],
  filename=c('handr_mekmut','handr_dnfgfr','foxf_ko','mesp_dnfgfr'),
  lfc=list(c(1,.5),c(1,.5),c(.75,.45),c(1,.7)),
  genes=list(
    prime.denovo[1:4],
    prime.denovo[1:4],
    prime.denovo[c(5,7,1:4)],
    prime.denovo[c(5,7,1:4)]
  ),
  col=list(
    c('red','blue','red','blue'),
    c('red','blue','red','blue'),
    c('forestgreen','gray28','red','blue','red','blue'),
    c('forestgreen','gray28','red','blue','red','blue')
  ),
  pch=list(
    c(19,19,1,1),
    c(19,19,1,1),
    c(19,19,19,19,1,1),
    c(19,19,19,19,1,1)
  ),
  ylim=list(
    c(-1.5,2.5),F,F, c(-2.5,2.6)
  ),
  MoreArgs=list(
    con=con,
    fdr=c(.05,.05),
    path='scatter',
    xlab='log2(RNAseq FC)',
    ylab='log2(ATACseq FC)'
  )
)
# Fig. 1E Spearman correlation
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
# Fig. 3A Spearman correlation
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


# Fig. 1F

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
# Fig. 1F Spearman correlation
cor(
  atac$condition_mesp_dnFGFR_vs_control[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
    ],
  atac$time_6hpf_vs_10hpf[
    union(peaksets$mespDep,peaksets$timeDep),'log2FoldChange'
  ],
  method = 'spearman'
)

# Fig. S10E
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
# Fig. S10E Spearman correlation
cor(
  atac$condition_mesp_dnFGFR_vs_control[
    peaksets$mespDep,'log2FoldChange'
    ],
  atac$condition_FoxF_KO_vs_control[
    peaksets$mespDep,'log2FoldChange'
  ],
  method = 'spearman'
)

