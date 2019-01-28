avgZ <- function(dev){
  sapply(
    unique(colData(dev)$condtime),
    function(x){
      z <- deviationScores(dev)[
        ,colData(dev)$condtime==x
      ]
      return(apply(z,1,mean,na.rm=T))
    }
  )
  # return(do.call(cbind,res))
}

writeChromVarHmap <- function(x,p=0.05,z=2){
  require(SummarizedExperiment)
  require(chromVAR)
  load(paste0(x,'/deviations.Rdata'))
  dev <- dev[,-6]
  hm <- chromVarHmap(dev,p,z)
  dir.eps(
    paste0(sub('.*\\/','',x),'fdr',as.character(p),'z',as.character(z)),
    height=nrow(hm@matrix)*.12+1,
    width=ncol(hm@matrix)*.25+3
  )
  draw(hm)
  dev.off()
}

chromVarHmap <- function(dev,p,z){
  require(ComplexHeatmap)
  require(circlize)
  require(TFBSTools)
  require(chromVAR)
  # load(paste0(x,'/deviations.Rdata'))
  sel <- !apply(deviations(dev),1,anyNA)
  family <- unlist(elementMetadata(dev)[["Family_Name"]])[sel]
  dev <- dev[sel,]
  
  sel <- !family%in%c(
    "ERF","AP2EREBP","BBRBPC","C2C2dof","Stat","ZFHD","Myb","MYB",
    "MYBrelated","ND","Trihelix","WRKY","POU,Homeobox",
    "promoter","NA"
  )
  dev <- dev[sel,]
  family <- family[sel]
  
  # dev <- dev[!apply(deviations(dev),1,anyNA),]
  colnames(dev) <- colData(dev)$Name
  row.names(dev) <- unlist(elementMetadata(dev)$DBID.1)
  
  diff_acc <- differentialDeviations(dev,'condtime')
  
  # sel <- order(diff_acc$p_value)[1:100]
  sel <- diff_acc$p_value_adjusted<p&!is.na(diff_acc$p_value_adjusted)
  dev <- dev[sel,]
  family <- family[sel]
  # dir.eps('top100',x,append.date = F,height=12)
  # draw(Heatmap(
  #   deviations(dev),
  #   cluster_columns = F,
  #   split=family,
  #   row_names_gp = gpar(cex=.5)
  # ))
  # dev.off()
  mat <- deviationScores(dev)
  mat[!is.finite(mat)] <- 0
  
  # average z-score
  mat <- do.call(
    cbind,
    lapply(
      split(
        row.names(colData(dev)),
        # define column order
        factor(
          colData(dev)$condtime,
          unique(colData(dev)$condtime)
        )
      ),
      function(x) apply(deviationScores(dev)[,x],1,mean,na.rm=T)
    )
  )
  
  sel <- apply(mat,1,function(x) any(abs(x)>z))
  # dir.eps('top100.z',x,append.date = F,height=12)
  return(Heatmap(
    mat[sel,],
    cluster_columns = F,
    split=family[sel],
    col=colorRamp2(c(-5,0,5),c('blue','white','red')),
    row_names_gp = gpar(cex=.5),
    row_title_rot = 0,
    row_title_gp = gpar(cex=.8),
    column_names_gp = gpar(cex=.8)
  ))
  dev.off()
  # return(split(elementMetadata(dev),family))
}

