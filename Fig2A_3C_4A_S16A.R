#Figs. 2A, 3C, 4A, S16A

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
  "MA_dnFGFR_LacZ_10hpf","FoxF10hpf_LacZ10hpf",
  "condtime_handrdnfgfr18hpf_handrlacz18hpf",
  "condtime_foxfcamras18hpf_handrlacz18hpf"
)]
atac <- getAtacLib(con,c(
  "condition_mesp_dnFGFR_vs_control","condition_FoxF_KO_vs_control",
  "condition_handr_dnFGFR_vs_control","condition_handr_MekMut_vs_control",
  "tissue_B7.5_vs_mesenchyme","condition_mesp_MekMut_vs_control"
))

prime.denovo <- scrna[c(
  'primedCardiac','primedASM','denovoCardiac','denovoASM','TVCP','STVC','ATM','mesenchyme'
)]
cols <- c(
  primedCardiac='red',primedASM='blue',
  denovoCardiac='hotpink',denovoASM='skyblue',
  TVCP='forestgreen',STVC='orange',ATM='gray28',
  ebfActivated='lightgreen',ebfInhibited='purple'
)

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
    scrna[13:16],
    scrna[13:16],
    scrna[c(7,10,13:16)],
    scrna[c(7,10,13:16)]
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
# Fig. 2A Spearman correlation
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
# Fig. 3C Spearman correlation
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

# Fig. S16A Spearman correlation
tmp <- mergeGenePeak(
  con,
  row.names(sig.sub(rna$condtime_foxfcamras18hpf_handrlacz18hpf)),
  row.names(sig.sub(atac$condition_handr_MekMut_vs_control,.5))
)
cor(
  rna$condtime_foxfcamras18hpf_handrlacz18hpf[tmp$GeneID,'log2FoldChange'],
  atac$condition_handr_MekMut_vs_control[tmp$PeakID,'log2FoldChange'],
  method='spearman'
)

