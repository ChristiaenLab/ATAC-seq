motifHmap <- function(
  con,hyper,dev,motifs,genes=NULL,file,p=.05,or=0
){
  require(ComplexHeatmap)
  require(circlize)
  require(TFBSTools)
  
  ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
  names(ma) <- as.character(seq(6,20,2))
  if(!is.null(genes)){
    ma <- ma[intersect(row.names(ma),genes),]
  }
  gene.names <- getGeneNames(con)
  # tf.family <- sapply(tags(motifs),'[[',"Family_Name")
  
  dat <- data.frame(
    TFGeneID=ID(motifs),
    TF=name(motifs),
    family=sapply(tags(motifs),'[[',"Family_Name"),
    row.names = names(motifs),
    stringsAsFactors=F
  )
  
  dat <- merge(dat,dev,'row.names')
  row.names(dat) <- dat$Row.names
  sel <- intersect(dat$Row.names,row.names(hyper[[1]]))
  dat <- dat[sel,-1]
  tf.kh.gene <- strsplit(dat$TFGeneID,';')
  dat <- dat[sapply(tf.kh.gene,function(x) any(x%in%genes)),]
  
  # motifs <- motifs[row.names(hyper[[1]])]
  
  odds <- sapply(hyper,'[',row.names(dat),'log2OddsRatio')
  fdr <- sapply(hyper,'[',row.names(dat),'p')
  fdr <- apply(fdr,2,p.adjust,'fdr')
  odds[fdr>p] <- 0
  colnames(odds) <- names(hyper)
  row.names(odds) <- row.names(dat)
  odds <- odds[
    !apply(odds<or,1,all),
    !apply(odds<or,2,all)
  ]
  
  dat <- merge(dat,odds,'row.names')[,-1]
  
  tf.kh.gene <- strsplit(dat$TFGeneID,';')
  times <- sapply(tf.kh.gene,length)
  dat <- mapply(
    function(x,y,z) cbind(
      GeneID=z,
      do.call(rbind,replicate(x,y,F)),
      stringsAsFactors=F
    ),
    times,
    split(dat,factor(dat$TFGeneID,dat$TFGeneID)),
    tf.kh.gene,
    SIMPLIFY = F
  )
  dat <- do.call(rbind,dat)
  
  dat <- merge(dat,ma,by.x="GeneID",by.y="row.names")
  
  # row.names(mat) <- make.unique(as.character(mat$TF))
  mat <- as.matrix(dat[,-1:-4])
  row.names(mat) <- khToName(dat$GeneID,ann$gene.names)
  
  heatmap_height <- nrow(mat)*.20+max(nchar(colnames(mat)))*.08
  lwidth <- max(nchar(dat$family))*.08
  rwidth <- max(nchar(row.names(mat)))*.08
  
  hm2 <- Heatmap(
    mat[,names(ma)],
    split=dat$family,
    heatmap_width = unit(ncol(ma)*.20+.5,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    cluster_columns = F,
    cluster_rows = T,
    cluster_row_slices = T,
    show_row_dend = F,
    row_dend_side='right',
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    column_gap=unit(.2,'in'),
    name="TF expression \n(log2FC vs. avg)"
  )
  
  hm1 <- Heatmap(
    mat[,colnames(odds)],
    split=dat$family,
    cluster_columns = F,
    cluster_row_slices=T,
    show_column_dend=F,
    show_row_dend=F,
    row_dend_side='right',
    heatmap_width = unit(rwidth+(ncol(odds))*.20+.5,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    col=colorRamp2(c(0,3),c('white','black')),
    column_gap=unit(.2,'in'),
    name='TF motif enrichment \n(log2OR)'
  )
  
  hm3 <- Heatmap(
    mat[,colnames(dev)],
    split=dat$family,
    cluster_columns = F,
    cluster_row_slices=T,
    show_column_dend=F,
    show_row_dend=F,
    row_dend_side='right',
    heatmap_width = unit(lwidth+(ncol(dev))*.20+.5,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    col = colorRamp2(c(-5,0,5),c("blue","white","red")),
    column_gap=unit(.2,'in'),
    name="TF motif accessibility \n(deviation z-score)"
  )
  
  dir.eps(
    file,
    width=(ncol(mat))*.20+4+lwidth+rwidth,
    height=heatmap_height+2
  )
  draw(hm3+hm2+hm1)
  dev.off()
}

