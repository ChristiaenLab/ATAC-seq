plotPeakMatches <- function(
  matches,file,split=NULL,
  col.split=NULL,
  ...
){
  require(ComplexHeatmap)
  require(circlize)
  matches <- t(matches)
  rowsel <- apply(matches,1,any)
  colsel <- apply(matches,2,any)
  split <- split[rowsel]
  col.split <- col.split[colsel]
  matches <- matches[rowsel,colsel,drop=F]
  
  if(ncol(matches)==0) return(NULL)
  
  motifs <- colnames(matches)
  matches <- t(apply(matches,1,as.numeric))
  colnames(matches) <- motifs
  height <- nrow(matches)/6+2
  width <- ncol(matches)*.25+5
  dir.eps(file,width=width,height=height)
  draw(col.hmap(
    matches,
    show_row_dend = F,
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    split=split,
    col.split=col.split,
    col = c('white','black'),
    ...
  ))
  dev.off()
}


peakHyper <- function(
  peaks,file,matches,split=NULL,
  fdr=T,p=0.05,
  col.split=NULL,
  # bg=Reduce(union,peaksets[c("timeDep","mespDep","handrDep")]),
  # width=6,height=8,
  ...
){
  require(ComplexHeatmap)
  require(circlize)
  peaksel <- sapply(
    peaks,
    function(x) row.names(matches)%in%x
  )
  row.names(peaksel) <- row.names(matches)
  hyper <- matchHyper(
    matches,
    peaksel
  )
  dir.tab(hyper,file)
  if(fdr) {
    fdr <- hyper[,grep("^fdr",colnames(hyper))]
  }else {
    fdr <- hyper[,grep("^p",colnames(hyper))]
  }
  rowsel <- apply(fdr<p,1,any)
  colsel <- apply(fdr<p,2,any)
  split <- split[rowsel]
  col.split <- col.split[colsel]
  fdr <- fdr[rowsel,colsel,drop=F]
  # fdr <- fdr[,colsel,drop=F]
  if(ncol(fdr)==0) return(NULL)
  hyper[!is.finite(hyper)] <- 0
  # fdrcol <- colorRamp2(
  #   c(0,-log2(0.05),-log2(.001)),c('dodgerblue4','white','firebrick4')
  # )
  height <- nrow(fdr)/6+2
  width <- ncol(fdr)*.25+5
  dir.eps(paste0(file,'FDR'),width=width,height=height)
  draw(col.hmap(
    -log2(fdr),
    # col=fdrcol,
    breaks=c(0,-log2(0.05),-log2(.001)),
    cluster_columns = F,show_row_dend = F,
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    split=split,
    col.split=col.split,
    ...
  ))
  dev.off()
  odds <- hyper[row.names(fdr),sub('fdr',"log2OddsRatio",colnames(fdr))]
  odds[fdr>p] <- 0
  # oddscol <- colorRamp2(
  #   seq(0,max(odds),length.out = 5),
  #   c("white","lightcyan",'steelblue1',"dodgerblue","navy")
  # )
  dir.eps(paste0(file,'Odds'),width=width,height=height)
  draw(abs.hmap(
    odds,
    # col=oddscol,
    cluster_columns = F,show_row_dend = F,
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    split=split,
    col.split=col.split,
    ...
  ))
  dev.off()
}

peakFamilyHyper <- function(
  peaks,file,familyID,matches,
  ...
){
  familyMatch <- sapply(
    unique(familyID),
    function(x) apply(
      matches[,familyID==x,drop=F],1,any
    )
  )
  peakHyper(peaks,file,familyMatch,...)
}
