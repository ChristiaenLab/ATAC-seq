motifHmap <- function(con,hyper,dev,motifs,file,p=.05,or=1.5){
  require(ComplexHeatmap)
  require(circlize)
  require(TFBSTools)
  
  odds <- sapply(hyper,'[',,'log2OddsRatio')
  fdr <- sapply(hyper,'[',,'padj')
  odds[fdr>p] <- 0
  motifs <- motifs[row.names(odds)]
  tf.family <- sapply(tags(motifs),'[[',"Family_Name")
  
  odds <- cbind(
    TFGeneID=ID(motifs),
    TF=names(motifs),
    family=tf.family,
    as.data.frame(odds),
    stringsAsFactors=F
  )
  odds <- odds[!apply(odds[,-1:-3]<or,1,all),c(rep(T,3),!apply(odds[,-1:-3]<or,2,all))]
  
  mat <- merge(odds,dev,'row.names')[,-1]
  
  tf.kh.gene <- strsplit(mat$TFGeneID,';')
  times <- sapply(tf.kh.gene,length)
  mat <- mapply(
    function(x,y,z) cbind(
      GeneID=z,
      do.call(rbind,replicate(x,y,F)),
      stringsAsFactors=F
    ),
    times,
    split(mat,factor(mat$TF,mat$TF)),
    tf.kh.gene,
    SIMPLIFY = F
  )
  # mat <- mapply(cbind,mat,GeneID=tf.kh.gene,stringsAsFactors=F,SIMPLIFY = F)
  # mat <- do.call(rbind,unlist(mat,F))
  # mat <- Reduce(function(x,y) {
  #   x[tf.kh.gene[[y]]] <- mat[[y]]
  #   return(x)
  # },1:length(mat),list())
  mat <- do.call(rbind,mat)
  
  
  ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
  names(ma) <- as.character(seq(6,20,2))
  mat <- merge(mat,ma,by.x="GeneID",by.y="row.names")
  row.names(mat) <- make.unique(as.character(mat$TF))
  # row.names(mat) <- gene.names[mat$GeneID,]
  
  # heatmap_width <- .25+max(nchar(row.names(mat)))*.05+max(nchar(as.character(mat$family)))*.05+(ncol(mat)-3)/4+heatmap_width
  heatmap_height <- nrow(mat)*.20+max(nchar(names(mat)))*.08
  lwidth <- max(nchar(mat$family))*.08
  rwidth <- max(nchar(row.names(mat)))*.08
  
  hm2 <- Heatmap(
    as.matrix(mat[,names(ma)]),
    split=mat$family,
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
    name="TF expression \n(log2FC)"
  )
  
  hm1 <- Heatmap(
    as.matrix(mat[,names(odds)[c(-1:-3)]]),
    # breaks = c(0,4),
    split=mat$family,
    cluster_columns = F,
    cluster_row_slices=T,
    show_column_dend=F,
    show_row_dend=F,
    row_dend_side='right',
    heatmap_width = unit(rwidth+(ncol(odds)-3)*.20+.5,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    # heatmap_width = unit((ncol(odds)-2)/4+heatmap_width,'in'),
    # heatmap_height = unit(heatmap_height,'in')
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    col=colorRamp2(c(0,3),c('white','black')),
    column_gap=unit(.2,'in'),
    name='TF motif enrichment \n(log2OR)'
  )
  
  hm3 <- Heatmap(
    as.matrix(mat[,colnames(dev)]),
    # breaks = c(0,4),
    split=mat$family,
    cluster_columns = F,
    cluster_row_slices=T,
    show_column_dend=F,
    show_row_dend=F,
    row_dend_side='right',
    heatmap_width = unit(lwidth+(ncol(dev))*.20+.5,'in'),
    heatmap_height = unit(heatmap_height,'in'),
    # heatmap_width = unit((ncol(odds)-2)/4+heatmap_width,'in'),
    # heatmap_height = unit(heatmap_height,'in')
    row_title_rot = 0,
    row_title_gp = gpar(cex=1),
    col = colorRamp2(c(-5,0,5),c("blue","white","red")),
    column_gap=unit(.2,'in'),
    name="TF motif accessibility \n(deviation z-score)"
  )
  
  dir.eps(
    file,
    width=(ncol(mat)-3)*.25+2+.25+lwidth+rwidth,
    height=heatmap_height+4
  )
  draw(hm3+hm2+hm1)
  dev.off()
}

