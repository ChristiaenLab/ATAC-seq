col.hmap <- function(
  # sets default colors for all heatmaps
  # adds option to split heatmap by column groups
  mat,...,
  breaks=NULL,
  center=median(mat),
  quant=c(.01,.99),
  col=c('blue','white','red'),
  col.split=NULL,
  show_row_names=T
){
  require(circlize)
  require(ComplexHeatmap)
  if(is.null(breaks)){
    breaks <- sapply(quant,function(x) quantile(mat,x))
  }
  if(length(breaks)<length(col)){
    if(!is.null(center)){
      breaks <- c(
        seq(breaks[1],center,length.out = (length(col)+1)/2),
        seq(center,breaks[2],length.out = (length(col)+1)/2)[-1]
      )
    }else{
      breaks <- seq(
        breaks[1],breaks[length(breaks)],length.out = length(col)
      )
    }
  }
  col <- colorRamp2(breaks,col)
  if(!is.null(col.split)) {
    mat <- lapply(split.data.frame(t(mat),col.split),t)
    # if heatmap is split by columns, only show row names for the last heatmap in the list
    if(show_row_names) show_row_names <- c(rep(F,length(mat)-1),T)
    hm <- mapply(
      Heatmap,
      mat,
      show_row_names=show_row_names,
      col_title=unique(col.split),
      MoreArgs=list(...,col=col)
    )
    hm <- Reduce('+',hm)
  } else hm <- Heatmap(
    mat,...,
    col=col,
    show_row_names=show_row_names
  )
  return(hm)
}

abs.hmap <- function(
  # default colors for heatmaps without negative values
  mat,...,
  breaks=c(0,max(mat)),
  quant=c(0,1),
  col=c("white","lightcyan",'steelblue1',"dodgerblue","navy"),
  center=NULL
){
  col.hmap(
    mat,...,
    breaks = breaks,
    quant=quant,
    col = col,
    center = center
  )
}

# center.hmap <- function(
#   mat,...,
#   breaks=c(
#     min(mat),max(mat)
#   ),
#   center=0,
#   quant=c(0,1),
#   col=c('dodgerblue4','white','firebrick4')
# ){
#   require(circlize)
#   require(ComplexHeatmap)
#   breaks <- c(seq(
#     breaks[1],breaks[2],length.out = length(col)
#   )
#   if(!is.null(quant)) {
#     quant <- seq(quant[1],quant[2],length.out = length(col))
#     breaks <- sapply(quant,function(x) quantile(mat,x,na.rm=T))
#   }
#   Heatmap(
#     mat,...,
#     col=colorRamp2(breaks,col)
#   )
# }

cor.heatmap <- function(
  # correlation heatmap for a matrix of read counts
  # columns should be replicates
  dat,cor.method='spearman',dist.method='spearman',
  cluster_rows=T,cluster_columns=T,
  dend_reorder=T,
  col=colorRamp2(c(
    quantile(atac.cor.sig,.05),
    quantile(atac.cor.sig,.25),
    quantile(atac.cor.sig,.75),
    quantile(atac.cor.sig,.95)
  ),c('dodgerblue4','steelblue','khaki1','firebrick')),
  moreArgsAnn=list(),...,cols=c('condition','time')
){
  # sel <- intersect(row.names(design),colnames(dat))
  # dat <- dat[,sel]
  # design <- design[sel,]
  atac.cor.sig <- apply(as.matrix(dat),2,function(x) apply(dat,2,cor,x,method=cor.method))
  # row.names(atac.cor.sig) <- design[row.names(atac.cor.sig),'Name']
  # colnames(atac.cor.sig) <- design[colnames(atac.cor.sig),'Name']
  ht <- col.hmap(
    atac.cor.sig,
    cluster_rows=cluster_rows,cluster_columns = cluster_columns,
    clustering_distance_columns = dist.method,
    clustering_distance_rows = dist.method,
    col=col,
    row_dend_reorder = dend_reorder,
    column_dend_reorder = dend_reorder,
    ...
  )
  # ann <- do.call(rowAnnotation,append(moreArgsAnn,list(design[,cols]),0))
  # ht <- ht+ann
  return(ht)
}

# accepts a counts matrix and a vector of element lengths
# computes RPKM for the matrix and passes the result to cor.heatmap
rpkm.heatmap <- function(dat,len,reads,...){
  cor.heatmap(
    dat/(len/1000)%*%(t(reads)/10^6),...
  )
}


# cor.heatmap <- function(
#   dat,design,cor.method='spearman',dist.method='spearman',
#   cluster_rows=T,cluster_columns=T,
#   dend_reorder=nrow(design):1,
#   moreArgsAnn=list(),...,cols=c('condition','time')
# ){
#   require(ComplexHeatmap)
#   require(circlize)
#   sel <- intersect(row.names(design),colnames(dat))
#   dat <- dat[,sel]
#   design <- design[sel,]
#   atac.cor.sig <- apply(dat,2,function(x) apply(dat,2,cor,x,method=cor.method))
#   row.names(atac.cor.sig) <- design[row.names(atac.cor.sig),'Name']
#   colnames(atac.cor.sig) <- design[colnames(atac.cor.sig),'Name']
#   ht <- Heatmap(
#     atac.cor.sig,
#     col=colorRamp2(c(
#       quantile(atac.cor.sig,.05),
#       quantile(atac.cor.sig,.95)
#     ),c('white','blue')),
#     cluster_rows=cluster_rows,cluster_columns = cluster_columns,
#     clustering_distance_columns = dist.method,
#     clustering_distance_rows = dist.method,
#     row_dend_reorder = dend_reorder,
#     column_dend_reorder = dend_reorder,
#     ...
#   )
#   ann <- do.call(rowAnnotation,append(moreArgsAnn,list(design[,cols]),0))
#   ht <- ht+ann
#   return(ht)
# }

rnaAtacCor <- function(con,rna,atac,ann,lfc=c(0,0),fdr=c(.05,.05),method="spearman"){
  # correlates RNA-seq and ATAC-seq lbraries
  # accepts a connection to a peak annotation database, a list of RNA-seq results from DESeq2,
  # a list of ATAC-seq results from DESeq2, a 2-element vector of log2FC cutoffs, 
  # a 2-element vector of FDR cutoffs, and a correlation method.
  # the first element of lfc and fdr are RNA-seq cutoffs
  # the second elements are ATAC-seq cutoffs
  # merges significant genes with significant peaks based on gene-to-peak annotation from con
  # returns correlation of log2FC values
  
  rna <- lapply(rna,sig.sub,lfc[1],fdr[1])
  atac <- lapply(atac,sig.sub,lfc[2],fdr[2])
  get.cor <- function(rna,atac) {
    sel <- mergeGenePeak(con,row.names(rna),row.names(atac))
    return(cor(
      rna[sel$GeneID,'log2FoldChange'],
      atac[sel$PeakID,"log2FoldChange"],
      method=method
    ))
  }
  cor <- sapply(rna,function(x) sapply(atac,function(y) get.cor(x,y)))
  return(cor)
}

# rnaAtacCor <- function(rna,atac,ann,lfc=c(0,0),fdr=c(.05,.05)){
#   # require(ComplexHeatmap)
#   rna <- lapply(rna,sig.sub,lfc[1],fdr[1])
#   atac <- lapply(atac,sig.sub,lfc[2],fdr[2])
#   rnasel <- unique(unlist(lapply(rna,row.names)))
#   atacsel <- unique(unlist(lapply(atac,row.names)))
#   peaks <- geneToPeak(
#     rnasel,
#     ann[atacsel,]
#   )
#   # peaks <- lapply(peaks,function(y) atac[y,,drop=F])
#   rna <- lapply(rna,'[',names(peaks),,drop=F)
#   rna <- lapply(rna,na.omit)
#   get.cor <- function(genes,atac) {
#     atac <- Filter(
#       function(x) length(x)>0,
#       lapply(
#         peaks[row.names(genes)],
#         function(z) na.omit(atac[z,'log2FoldChange'])
#       )
#     )
#     if(length(atac)==0) return(0)
#     genes <- genes[names(atac),'log2FoldChange']
#     tmp <- do.call(rbind, mapply(
#       cbind, genes,atac
#     ))
#     return(cor(tmp[,1],tmp[,2],method='spearman'))
#   }
#   cor <- sapply(rna,function(x) sapply(atac,function(y) get.cor(x,y)))
#   return(cor)
# }

de.daHeatmap <- function(
  rna,atac,genes,ann,gene.names,levels=NULL,col=NULL,p=c(.05,.05),lfc=c(1,.5),...
){
  rna <- lapply(rna,sig.sub,lfc[1],p[1])
  atac <- lapply(atac,sig.sub,lfc[1],p[1])
  combn()
  require(ComplexHeatmap)
  dat <- mergeRnaAtac(list(rna),list(atac),genes,ann,gene.names,p,lfc)
  
  if(!is.null(levels)&!is.null(col)) {
    levels <- as.data.frame(do.call(rbind,mapply(cbind,levels,names(levels))),stringsAsFactors=F)
    levels <- levels[!duplicated(levels[,1]),,drop=F]
    row.names(levels) <- levels[,1]
    levels$col <- col[levels[,2]]
  }
  
  atac <- do.call(rbind,mapply(
    function(x,y) cbind(x,row.names(y),y),
    x=names(dat$peaks),
    y=dat$peaks
  ))
  atac <- cbind(atac,dat$rnasel[atac[,1],c('gene','split')])
  atac <- as.matrix(atac)
  sapply(
    atac[duplicated(atac[,2]),2],
    function(x) {
      sel <- atac[,2]==x
      atac[sel,1] <<- do.call(paste,append(atac[sel,1],list(sep=';')))
      atac[sel,'gene'] <<- do.call(paste,append(atac[sel,'gene'],list(sep=';')))
      atac[sel,'split'] <<- do.call(paste,append(unique(atac[sel,'split']),list(sep=';')))
    }
  )
  atac <- atac[!duplicated(atac),]
  atac <- cbind(atac,levels[atac[,1],2])
  # atac[is.na(atac[,6]),6] <- 'none'
  # atac[is.na(atac[,7]),7] <- 'white'
  counts <- as.matrix(counts[atac[,2],])
  row.names(counts) <- atac[,4]
  # row.sd <- apply(counts,1,sd)*sqrt((ncol(counts)-1)/ncol(counts))
  # row.mean <- apply(counts,1,mean)
  # row.z <- (counts-row.mean)/row.sd
  ht <- Heatmap(
    counts,
    # row.z,
    split = atac[,'split'],
    # row_names_gp = gpar(cex=cex.row),
    cluster_columns = F
    # show_row_names = F,...
  )
  ann <- rowAnnotation(data.frame(x=atac[,6]),'scRNA',list(x=col))
  hmap <- ht+ann
  return(hmap)
}
