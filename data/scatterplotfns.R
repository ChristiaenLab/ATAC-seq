scatterRnaAtac <- function(
  con,rna,atac,
  genes=NULL,col=NULL,pch=NULL,
  filename,path='.',
  lfc=c(1,1),fdr=c(0.05,0.05),
  cex=1.5,
  xlim=F,ylim=F,
  ...
){
  de <- geneToPeak(con,row.names(sig.sub(rna,lfc[1],fdr[1])))
  da <- mergeGenePeak(con,de$GeneID,row.names(sig.sub(atac,lfc[2],fdr[2])))
  
  dat <- data.frame(
    rna[de$GeneID,'log2FoldChange'],
    atac[de$PeakID,'log2FoldChange']
  )
  sig <- data.frame(
    rna[da$GeneID,'log2FoldChange'],
    atac[da$PeakID,'log2FoldChange']
  )
  da.genesets <- sapply(genes,intersect,da$GeneID)
  da.pts <- lapply(da.genesets,function(x) sig[da$GeneID%in%x,])
  
  if(length(xlim)<2) {xlim <- range(dat[,1],na.rm = T)}
  if(length(ylim)<2) {ylim <- range(sig[,2],na.rm = T)}
  
  dir.eps(filename,path)
  plot(NULL,xlim=xlim,ylim=ylim,cex.axis=cex,...)
  points(dat,col='gray',cex=cex,pch=19)
  points(sig,cex=cex,pch=1)
  mapply(points,da.pts,col=col,cex=cex,pch=pch)
  segments(xlim[1],0,xlim[2],0)
  segments(0,ylim[1],0,ylim[2])
  
  gene.names <- getGeneNames(con)[da$GeneID,"UniqueNAME"]
  gene.names <- sub(".*_",'',gene.names)
  gene.names <- sub("KH2013:","",gene.names)
  
  text(sig,cex=.5,labels = gene.names)
  
  dev.off()
}

