library(DBI)
source('data/fgsea.R')
source('data/DESeqFns.R')
source('data/sqlfns.R')
source('data/dirfns.R')

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

# write GSEA results to database
gsea <- mapply(
  get.fgsea,
  atac,
  names(atac),
  MoreArgs=list(
    peaklists=peaks,
    p=1,path='gsea'
  ),
  SIMPLIFY = F
)

dbWriteTable(con,'gsea',melt.rename(gsea,'comparsion'))
