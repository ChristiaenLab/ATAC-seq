library(DBI)
library(DESeq2)
source('fgsea.R')
source('data/sqlfns.R')
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
peakGeneAnnotation <- getAnnotation(con)
bulkGS <- getBulkRNA(con)
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
atac <- getAtac(con)
feat <- apply(
  peakGeneAnnotation$features[,-length(peakGeneAnnotation$features)],
  2,
  function(x) row.names(peakGeneAnnotation$features)[x]
)
feat$promoter <- union(feat$promoter1k,feat$promoter500)

peaks <- append(lapply(
  Reduce(append,list(scrna,bulkGS,feat)),
  function(x)unique(unlist(peakGeneAnnotation$geneToPeak[x]))
),peaksets)

path <- 'gsea'
mapply(get.fgsea,atac,names(atac),MoreArgs=list(peaklists=peakGeneAnnotation$geneToPeak,p=1,path=paste0(path,'/genes')))
mapply(get.fgsea,atac,names(atac),MoreArgs=list(peaklists=peakGeneAnnotation$geneToPeak,path=paste0(path,'/genesSig')))

mapply(get.fgsea,atac,names(atac),MoreArgs=list(
  peaklists=peaks,
  p=1,path=paste0(path,'/peaksets')
))
mapply(get.fgsea,atac,names(atac),MoreArgs=list(
  peaklists=lapply(scrna,function(x)unique(unlist(peakGeneAnnotation$geneToPeak[x]))),
  path=paste0(path,'/peaksetsSig')
))

