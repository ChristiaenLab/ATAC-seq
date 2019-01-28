get.dds <- function(
  dat,expDesign,design,path='.',comparisons=NULL,cooksCutoff=T,parallel=T,alpha=.1,independentFiltering=T,...
){
  require(DESeq2)
  expDesign[expDesign==T] <- "yes"
  expDesign[expDesign==F] <- "no"
  
  filt <- dat[,row.names(expDesign)]
  dds <- DESeqDataSetFromMatrix(
    countData=filt, 
    colData = expDesign, 
    design = design)
  dds <- DESeq(dds,...,parallel = parallel)
  # write normalized counts
  # this should be identical for either DESeq test
  normcts <- as.data.frame(counts(dds,normalized=T))
  dir.csv(normcts,'counts_norm',path)
  
  save(dds,file = mkdate('dds','Rdata',path))
  if(!is.null(comparisons)){
    return(res.list(
      dds,comparisons,path,
      cooksCutoff=cooksCutoff,parallel=parallel,alpha=alpha,
      independentFiltering=independentFiltering
    ))
  }
}

res.list <- function(
  dds,comparisons,path='.',format = "DataFrame",alpha = .1,
  cooksCutoff=T,parallel=T,independentFiltering=T
){
  require(DESeq2)
  res <-  mapply(
    results,list(dds),comparisons,
    format = format,alpha = alpha,parallel=parallel,
    cooksCutoff=cooksCutoff,
    independentFiltering=independentFiltering
  )
  names(res) <- sapply(res,function(x) gsub(
    '\\+','',gsub('\\s','_',sub('.*:\\s','',elementMetadata(x)$description[2]))
  ))
  save(res,file = mkdate('res','Rdata',path))
  return(res)
}


is.sig <- function(x,lfc=1,p=0.05,tail='both',na.as.false=T) {
  res <- x$padj<p
  if(length(lfc)==2){
    res <- (x$log2FoldChange<lfc[1]|x$log2FoldChange>lfc[2])&res
  }else if(tail=="both"){
    res <- abs(x$log2FoldChange)>lfc&res
  }else if(tail=='upper'){
    res <- x$log2FoldChange>lfc&res
  }else if(tail=='lower'){
    res <- x$log2FoldChange< -lfc&res
  }else if(lfc<0){
    res <- x$log2FoldChange<lfc&res
  }else res <- x$log2FoldChange>lfc&res
  if(na.as.false) res <- res&!is.na(res)
  return(res)
}

which.sig <- function(x,lfc=1,p=0.05,tail='both') which(is.sig(x,lfc,p,tail))
sig.sub <- function(x,lfc=1,p=0.05,tail='both') x[which(is.sig(x,lfc,p,tail)),]
any.sig <- function(x,lfc=1,p=0.05,tail='both') if((length(lfc)==length(x)&is.character(tail))|is.list(lfc)){
  apply(mapply(is.sig,x,lfc,p,tail),1,any,na.rm=T)
}else {
  apply(sapply(x,is.sig,lfc,p,tail),1,any,na.rm=T)
}