library(DBI)
library(ComplexHeatmap)
source("data/dirfns.R")
source("data/sqlfns.R")
source("data/cor.heatmap.R")
source("data/DESeqFns.R")

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
peaksets <- getPeaksets(con)

len <- dbGetQuery(con,'SELECT PeakID,end-start+1 FROM peakfeature',row.names="PeakID")
dat <- dbReadTable(con,'atacreads',row.names="PeakID")
expDesign <- dbGetQuery(
  con,
  'SELECT * FROM ataclib WHERE condition="control" AND (NOT omit OR Name="gfp_10hpf_2") AND NOT time IN ("15hpf","20hpf") AND NOT KO AND tissue!="Ef1a" AND NOT Name IN ("gfp_10hpf_3","LacZ_18hpf_2","LacZ_10hpf_3")',
  row.names='lib'
)
dat <- dat[ ,row.names(expDesign) ]
names(dat) <- expDesign$Name

dir.eps('rpkmSpearmanReorder')
draw(rpkm.heatmap(
  dat,unlist(len),expDesign$reads,
  dend_reorder = c(rep(1,6),rep(3,5),rep(2,3)),
  col=colorRamp2(seq(0.1,.9,length.out = 5),c("white","lightcyan",'steelblue1',"dodgerblue","navy"))
))
dev.off()

dat <- dat[
  do.call(union,setNames(peaksets[c('timeDep','tissueDep')],NULL)),
]
len <- len[do.call(union,setNames(peaksets[c('timeDep','tissueDep')],NULL)),]

dir.eps('rpkmSpearmanReorderSig')
draw(rpkm.heatmap(
  dat,len,expDesign$reads,
  dend_reorder = c(rep(1,6),rep(3,5),rep(2,3)),
  col=colorRamp2(seq(0,1,length.out = 5),c("white","lightcyan",'steelblue1',"dodgerblue","navy"))
))
dev.off()

rna <- getRnaDat(con)[-1]
atac <- getAtacLib(con,c(
  "condition_handr_MekMut_vs_control",
  "condition_handr_dnFGFR_vs_control",
  "condition_mesp_dnFGFR_vs_control",
  "condition_FoxF_KO_vs_control",
  "time_6hpf_vs_10hpf"
))

dir.eps("rnaAtacCor")
draw(col.hmap(
  rnaAtacCor(con,rna,atac),cluster_rows = F,cluster_columns = F,center=0
  # col=colorRamp2(seq(-0.5,.5,length.out = 5),c("white","lightcyan",'steelblue1',"dodgerblue","navy"))
))
dev.off()
