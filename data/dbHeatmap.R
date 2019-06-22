dbHeatmap <- function(
  con,rna,atac,gene.sets,
  peak.sets,peak.cols,
  more.genes,more.gene.cols,
  file=NULL,path='.',
  peak.pch=18,
  peak.cex=1.2,
  gene.peak.intersect=T,
  p=c(0.05,0.05),
  cex.axis=1.2,...
) {
  require(DBI)
  require(grid)
  require(gridExtra)
  require(ComplexHeatmap)
  require(circlize)
  require(reshape2)
  
  gene.names <- dbReadTable(con,'gene_name',row.names="GeneID")
  if(gene.peak.intersect){
    gene.sets <- lapply(
      gene.sets,
      function(x) unique(mergeGenePeak(
        con,x,unique(unlist(peak.sets))
      )$GeneID)
    )
  }
  
  gene.peak <- lapply(
    lapply(gene.sets,geneToPeak,con=con),
    do.call,
    what=function(GeneID,PeakID) split(PeakID,GeneID)
  )
  
  dat <- as.data.frame(do.call(rbind,mapply(
    cbind,
    lapply(gene.peak,names),
    names(gene.peak),
    SIMPLIFY = F
  )),stringsAsFactors=F)
  names(dat) <- c("GeneID","geneset")
  dat <- dat[!duplicated(dat$GeneID),]
  row.names(dat) <- dat$GeneID
  dat$gene <- gene.names[dat$GeneID,"UniqueNAME"]
  dat$gene <- sub('KH2013:','',dat$gene)
  dat$gene <- sub('.*_','',dat$gene)
  dat$gene <- sub('(.{,12}).*','\\1',dat$gene)
  dat$gene <- sapply(dat$gene,function(x) if(
    grepl("^[A-Z]{2}",x)&!(grepl("^KH\\.",x)|grepl("^SI:",x))
  ){
    tmp <- sub("(^.)([A-Z]+)",'\\2',x)
    return(sub(tmp,tolower(tmp),x))
  } else {x})
  
  ngenes <-  unlist(table(dat$geneset))
  
  rna <- sapply(rna,'[',dat$GeneID,'log2FoldChange')
  row.names(rna) <- dat$GeneID
  
  ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
  names(ma) <- as.character(seq(6,20,2))
  
  comb <- cbind(ma[dat$GeneID,],rna[dat$GeneID,,drop=F])
  comb[is.na(comb)] <- 0
  row.names(comb) <- dat$GeneID
  # sel <- row.names(dat$rna)[!row.names(dat$rna)%in%row.names(ma)]
  # ma[setdiff(unlist(gene.sets),row.names(ma)),] <- 0
  
  # ma <- ma[dat$GeneID,]

  ann.cols <- melt(more.genes)
  ann.cols <- split(ann.cols,ann.cols[,3])
  ann.cols <- lapply(ann.cols,function(x){
    x <- x[-duplicated(x[,1]),]
    fill <- setdiff(row.names(comb),x[,1])
    x <- as.matrix(x)
    if(length(fill)>0){
      x <- rbind(x, cbind(fill,NA,NA)
      )
    }
    row.names(x) <- x[,1]
    x <- x[row.names(comb),]
    return(x)
  })
  ann.cols <- as.data.frame(sapply(ann.cols,'[',,2))
  ann.cols <- mapply(function(x,y){
    x[!x%in%names(y)] <- NA
    return(x)
  },ann.cols,more.gene.cols)
  gene.ann <- rowAnnotation(
    df=as.data.frame(ann.cols),
    col=more.gene.cols,
    na_col='white'
  )
  name.ann <- rowAnnotation(
    text=row_anno_text(
      dat$gene
    )
  )

  heatmap_width <- ncol(comb)/4+length(atac)*2.5+.25+ncol(ann.cols)/4
  heatmap_height <- nrow(comb)/4+max(nchar(c(names(comb),names(gene.sets))))*.05+.5

  peakdat <- geneToPeak(con,dat$GeneID)
  peakdat <- split(peakdat$PeakID,peakdat$GeneID)
  atac <- lapply(atac, function(x) lapply(peakdat,function(y) x[y,]))
  # atac <- lapply(atac,'[',peakdat$PeakID,)
  # atac <- lapply(atac,split,peakdat$GeneID)
  atac <- lapply(atac,'[',row.names(comb))
  
  rightAnn <- lapply(
    atac,
    atac.ann,
    comb,
    peak.sets,
    peak.cols,
    peak.pch,
    peak.cex
  )
  # names(rightAnn) <- names(atac)
  rightAnn <- do.call(HeatmapAnnotation,append(rightAnn,list(which='row')))
  # heatmap_width2 <- length(rna)/5
  
  hm <- Heatmap(
    as.matrix(comb),
    split = dat$geneset,
    cluster_columns = F,
    show_row_names = F,
    gap=unit(4,'mm'),
    column_gap=unit(4,'mm'),
    heatmap_width = unit(heatmap_width,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    column_split = factor(c(
      rep('time',ncol(ma)),
      rep('condition',ncol(rna))
    ),c('time','condition')),
    right_annotation = rightAnn,
    cluster_rows = T,
    cluster_row_slices = F,
    col = colorRamp2(c(-4,0,4),c('blue','white','red'))
  )
  
#   hm2 <- Heatmap(
#     sapply(rna,'[',dat$GeneID,'log2FoldChange'),
#     cluster_columns = F,
#     cluster_rows = F,
#     na_col = 'white',
#     right_annotation = rightAnn,
#     heatmap_width = unit(heatmap_width2+length(atac)*2.5,'in'),
#     heatmap_height = unit(heatmap_height,'in')
#   )
  res <- hm+gene.ann+name.ann
  width <- heatmap_width+4
  height <- heatmap_height+2
  if(is.character(file)){
    dir.eps(file,path,width=width,height=height)
    draw(res)
    dev.off()
  }
}

filtpts <- function(x,y,sel,col='gray',cex=1.2,pch=18,...){
  y <- y[sel]
  x <- x[sel]
  if(length(y)>0) return(grid.points(
    x,y,pch=pch,
    gp=gpar(col=col,cex=cex,...), 
    default.units = "native"
  ))
}
    
atac.ann.fn <- function(index, k, n) {
  
  x <- length(index):1-.5
  # atac.sub <- subset_vector(atac,index)
  # atac.names <- atac.names[index]
  atac.sub <- atac[index]
  atac.y <- lapply(atac.sub,'[',,'log2FoldChange')
  atac.x <- mapply(
    function(z,w) jitter(rep(z,length(atac.y[[w]])),amount=.23),
    x,1:length(index)
  )
  
  yscale=c(-.05,length(atac.sub)+.05)
  pushViewport(viewport(
    yscale=yscale,
    xscale = range(unlist(lapply(atac,'[',,'log2FoldChange')))+c(-.25,.25)
  ))
  grid.rect()
  sapply(x,grid.abline,0)
  mapply(
    filtpts,
    y=atac.x,
    x=atac.y,
    sel=T
  )
  peak.sel <- lapply(
    peak.sets,
    function(i) lapply(atac.sub, function(j) row.names(j)%in%i)
  )
  mapply(
    function(
      peak.sel,peak.cols,peak.pch,peak.cex
    ) mapply(
      filtpts,
      x=atac.y,
      y=atac.x,
      sel=peak.sel,
      MoreArgs = list(
        col=peak.cols,pch=peak.pch,cex=peak.cex
      )
    ),
    peak.sel,
    peak.cols,
    peak.pch,
    peak.cex
  )
  # grid.text(gene.names[names(atac.sub),'UniqueNAME'],0,x,default.units = 'native')
  # grid.text(as.character(index),0,x,default.units = 'native')
  grid.lines(
    0,yscale,
    default.units = "native",
    gp=gpar(lty = 2)
  )
  if(k == n) grid.xaxis(draw=T)
  
  popViewport()
}

atac.ann <- function(
  atac,ma,peak.sets,peak.cols,peak.pch,peak.cex
) AnnotationFunction(
  fun = atac.ann.fn,
  var_import = list(
    # atac.names=names(atac),
    atac=atac,
    filtpts=filtpts,
    peak.sets=peak.sets,
    peak.cols=peak.cols,
    peak.pch=peak.pch,
    peak.cex=peak.cex,
    gene.names=gene.names
  ),
  n = nrow(ma),
  subsetable = TRUE,
  # subset_rule = subset_vector,
  width = unit(2.5, "in"),
  # height = unit(2, "cm"),
  which='row'
)

