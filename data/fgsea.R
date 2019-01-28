dir.fgsea <- function(
  deseq='2018-05-08/peakome/DA_LRT',peaklists='2018-03-13/peakome/scRNApeaks'
) {
  peaklists <- lrtab(peaklists,read.bed,'\\.bed')
  load(paste0(deseq,'/res.Rdata'))
  lapply(
    names(res),
    function(x) get.fgsea(res[[x]],x,peaklists)
  )
}

get.fgsea <- function(atac,file,peaklists,p=1,lfc=0,minSize=5,path='gsea',plotset=F,nperm=1000){
  require('fgsea')
  require('ggplot2')
  
  path <- paste0(path,'/',file)
  
  atac <- get.sig(atac,lfc=lfc,p=p)
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
      nperm=1000
    ), ranks=ranks,peaklists=peaklists
  )
  save(fgseaRes,file = mkdate('fgseaRes','Rdata',path))
  
 # if(length(fgseaRes$peaklists)<100&length(fgsea$peaklists)>0&plotset){
 #   lapply(names(fgseaRes$peaklists),function(x) dir.gg(
 #     plotEnrichment(fgseaRes$peaklists[[x]], fgseaRes$ranks),x,path
 #   ))
 #   dir.pdf('gseaTable',path)
 #   plotGseaTable(fgseaRes$peaklists,fgseaRes$ranks,fgseaRes$fgsea,gseaParam = 0.5)
 #   dev.off()
 # }
}
