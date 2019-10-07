get.dds <- function(
  dat,expDesign,design,path='.',comparisons=NULL,cooksCutoff=T,
  parallel=T,alpha=.1,independentFiltering=T,...
){
  # wrapper function for DESeq2 which saves the dataset and results object to a specified directory
  # comparisons should be a list of comparisons passed to res.list
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
  res <- NULL
  if(!is.null(comparisons)){
    res <- res.list(
      dds,comparisons,path,
      cooksCutoff=cooksCutoff,parallel=parallel,alpha=alpha,
      independentFiltering=independentFiltering
    )
  }
  return(list(dds=dds,res=res))
}

res.list <- function(
  dds,comparisons,path='.',format = "DataFrame",alpha = .1,
  cooksCutoff=T,parallel=T,independentFiltering=T
){
  # accepts a DESeqDataSet and a list of comparisons
  # applies results to the DESeqDataSet for each comparison in comparisons
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
  # returns a logical vector of which rows in x are significant
  # x     a data.frame with columns p.adjust and log2FoldChange
  # lfc   log2FoldChange cutoff 
  #       If lfc is a vector of length 2, the first element is taken as the lower limit 
  #       and the second element is taken as the upper limit.
  # p     FDR cutoff
  # tail  may be one of "both","upper",or "lower".
  #       If "both", returns TRUE if x > lfc or x < -lfc
  #       If "upper",returns TRUE if x > lfc
  #       If "lower", returns TRUE if x < -lfc
  # na.as.false logical indicating whether to return FALSE if is.na(x)
  res <- x$padj<=p
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

# wrapper function for converting is.sig to numeric indices
which.sig <- function(x,lfc=1,p=0.05,tail='both') which(is.sig(x,lfc,p,tail))

# wrapper function for subsetting x by is.sig(x)
sig.sub <- function(x,lfc=1,p=0.05,tail='both') x[which(is.sig(x,lfc,p,tail)),]

# applies is.sig to a list of data.frames
# returns the union of all significant rows in any data.frame in the list
# lfc may be a list or a vector of the same length as x, 
# in which case each element in lfc is taken as a separate cutoff
any.sig <- function(x,lfc=1,p=0.05,tail='both') if((length(lfc)==length(x)&is.character(tail))|is.list(lfc)){
  apply(mapply(is.sig,x,lfc,p,tail),1,any,na.rm=T)
}else {
  apply(sapply(x,is.sig,lfc,p,tail),1,any,na.rm=T)
}