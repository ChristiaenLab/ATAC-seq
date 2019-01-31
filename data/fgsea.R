get.fgsea <- function(atac,file,peaklists,p=1,lfc=0,minSize=5,path='gsea',nperm=10000){
  # wrapper function for fgsea
  # atac        a data.frame of differential accessibility results with columns "log2FoldChange" and "padj"
  # file        output subdirectory to save results
  # peaklists   a list of PeakID vectors treated as gene sets
  # p           padj cutoff for filtering atac
  # lfc         log2FoldChange cutoff for filtering atac
  # minSize     minimum length for filtering vectors in peaklists
  # path        top-level directory to write file
  # nperm       Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
  require('fgsea')
  require('ggplot2')
  
  path <- paste0(path,'/',file)
  
  atac <- sig.sub(atac,lfc=lfc,p=p)
  ranks <- atac$log2FoldChange
  names(ranks) <- row.names(atac)
  peaklists <- sapply(peaklists,function(x) x[x%in%names(ranks)])
  peaklists <- Filter(function(x) length(x)>minSize,peaklists)
  
  fgseaRes <- list(
    fgsea=fgsea(
      pathways = peaklists, 
      stats = ranks[order(ranks)],
      minSize=minSize,
      maxSize=length(ranks)-1,
      nperm=nperm
    ), 
    ranks=ranks,
    peaklists=peaklists
  )
  save(fgseaRes,file = mkdate('fgseaRes','Rdata',path))
  return(as.data.frame(fgseaRes$fgsea)[,-ncol(fgseaRes$fgsea)])
}
