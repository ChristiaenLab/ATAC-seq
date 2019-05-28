dbDiamond <- function(con,genes,peaks,file='diamond'){
  genes <- strsplit(genes,';')[[1]]
  genes <- lapply(genes,function(x) dbGetQuery(con,x)$GeneID)
  names(genes) <- as.character(1:length(genes))
  
  peaks <- strsplit(peaks,';')[[1]]
  peaks <- lapply(peaks,function(x) dbGetQuery(con,x)$PeakID)
  names(peaks) <- as.character(1:length(peaks))
  
  peakGeneAnnotation <- getAnnotation(con)
  gene.names <- dbReadTable(con,'gene_name',row.names="GeneID")
  
  masub <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
  names(masub) <- as.character(seq(6,20,2))
  
  # bulkGS <- dbReadTable(con,'bulkRNAgenesets',row.names='geneset')
  # bulkGS <- apply(bulkGS,1,function(x) dbGetQuery(con,x)$GeneID)
  
  # scrna <- sapply(c(
  #   'ASM','ebfActivated','ebfInhibited',
  #   'FHP','FHP14','Pancardiac',"SHP",
  #   "STVC","TVCP","ATM_genes_from_ANISEED"
  # ),function(x) dbReadTable(con,x,row.names="GeneID"))
  # 
  # prime.denovo <- list(
  #   primedCardiac=Reduce(union,sapply(
  #     scrna[c("Pancardiac","FHP","FHP14","SHP")],
  #     function(x) row.names(x)[x$Type=="Primed"]
  #   )),
  #   primedASM=row.names(subset(scrna$ASM,Type=="Primed")),
  #   denovoCardiac=Reduce(union,sapply(
  #     scrna[c("Pancardiac","FHP","FHP14","SHP")],
  #     function(x) row.names(x)[x$Type=="De Novo"]
  #   )),
  #   denovoASM=row.names(subset(scrna$ASM,Type=="De Novo")),
  #   TVCP=row.names(scrna$TVCP), 
  #   STVC=setdiff(row.names(scrna$STVC),row.names(scrna$TVCP)),
  #   ATM=setdiff(
  #     row.names(scrna$ATM_genes_from_ANISEED),
  #     Reduce(union,scrna[c("ASM","FHP","FHP14","SHP","Pancardiac","STVC","TVCP")])
  #   )
  # )
  # 
  # ebf <- sapply(scrna[c('ebfActivated','ebfInhibited')],row.names)
  
  atacdat <- mapply(dbReadTable,name=c(
    "condition_handr_dnFGFR_vs_control","condition_handr_MekMut_vs_control",
    "condition_mesp_dnFGFR_vs_control","condition_FoxF_KO_vs_control","tissue_B7_5_vs_mesenchyme"
  ),MoreArgs = list(conn=con,row.names="PeakID"),SIMPLIFY = F)
  
  rnadat <- mapply(dbReadTable,name=c(
    "condtime_handrdnfgfr18hpf_handrlacz18hpf",
    "condtime_foxfcamras18hpf_handrlacz18hpf",
    "MA_dnFGFR_LacZ_10hpf","FoxF10hpf_LacZ10hpf"
  ),MoreArgs = list(conn=con,row.names="GeneID"),SIMPLIFY = F)
  rnadat$gfp.lacz <- rnadat$FoxF10hpf_LacZ10hpf
  rnadat$gfp.lacz$log2FoldChange <- 0
  rnadat$gfp.lacz$padj <- 1
  
  rnaAtacPanelHmap(
    genes= genes,
    file= file,
    which.sig=NULL,
    which.da.lib = function(x,y,z,w,v) rep(T,length(x)),
    which.da.gene = function(...) T,
    peak.sets= peaks,
    peak.cols= rep(2,length.out=length(peaks)),
    rna=rnadat,atac=atacdat,
    ma=masub,
    # levels = list(prime.denovo),
    ann=peakGeneAnnotation,
    gene.names=gene.names,
    # col=cols,
    lfc=list(
      c(1,1,1,.75,1),
      rep(.5,5)
    )
  )
}

dbMergeGenePeak <- function(con,rna,atac,gene.sets,peak.sets){
  gene.peak <- dbMergeGenePeak
}

# con is a database
# rna is a vector of rnaseq result table names in con
# atac is a vector of atacseq result tables corresponding to rna
# gene.sets is a vecor of SQL queries returning the GeneIDs used to split the plot panels
# peak.sets is a vector of SQL queries returning the PeakIDs used to color code the diamonds
# more.genes is a list of vectors of SQL queries returning GeneIDs. 
# Each list will be used to add a gene color key below the panels 
# if gene.peak.intersect is TRUE, only genes with a peak in any of the peak.sets will be plotted.
# 

dbDiamondplot <- function(
  con,rna,atac,gene.sets,
  peak.sets,peak.cols,
  more.genes,more.gene.cols,
  file=NULL,path='.',gene.peak.intersect=T,p=c(0.05,0.05),cex.axis=1.2,...
) {
  require(DBI)
  require(grid)
  require(gridExtra)
  
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
  
  
  ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
  names(ma) <- as.character(seq(6,20,2))
  # sel <- row.names(dat$rna)[!row.names(dat$rna)%in%row.names(ma)]
  ma[setdiff(unlist(gene.sets),row.names(ma)),] <- 0
  dend.order <- lapply(
    gene.peak,
    function(x){
      ma <- ma[names(x),]
      dend <- as.dendrogram(hclust(dist(as.matrix(ma))))
      return(row.names(ma)[order.dendrogram(dend)])
    }
  )
  
  gene.peak <- mapply('[',gene.peak,dend.order,SIMPLIFY = F)
  
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
  
  ngenes <-  unlist(table(dat$geneset))
  
  ma <- ma[dat$GeneID,]
  
  # rnaplot <- mapply(
  #   diamondPlotPanel,
  #   rna,atac,
  #   MoreArgs = list(gene.peak,peak.sets,peak.cols)
  # )
  # atacplot <- mapply(
  #   diamondPlotPanel,
  #   rna,atac,
  #   MoreArgs = list(
  #     gene.peak,peak.sets,peak.cols,
  #     fn=panel.atac,ylim='atac'
  #   )
  # )
  # plots <- mapply(doubleYScale,atacplot,rnaplot,SIMPLIFY = F)
  plots <- mapply(
    diamondPlotPanel,
    rna,atac,
    MoreArgs = list(
      gene.peak,peak.sets,peak.cols
    )
  )
  # return(plots)
  
  zmax <- quantile(abs(unlist(ma)),.99)
  ma[ma>zmax] <- zmax
  ma[ma< -zmax] <- -zmax
  at <- seq(-zmax,zmax,length.out = 16)
  lplots <- levelplot(
    z~x*y|geneset,
    data = data.frame(
      expand.grid(
        x=1:nrow(ma),#factor(row.names(ma),row.names(ma)),
        y=names(ma)
      ),
      z=unlist(ma),
      geneset=dat$geneset
    ),
    layout=c(length(gene.sets),1),
    scales=list(
      x=list(
        rot=90,relation='free',cex=0,draw=F,axs='i'
      ),
      y=list(
        cex=cex.axis,draw=T,relation='free'
      )
    ),
    between=list(x=2),                       #between creates more space
    # between the panels - you may need to adjust this value
    par.settings=list(layout.widths=list(right.padding=6)), #this creates
    ylab=NULL,xlab=NULL,
    prepanel=function(x, y, subscripts,...) {
      xsub <- unique(x[subscripts])
      return(list( 
        xat = as.numeric(xsub)
      ))
    },
    col.regions = colorRampPalette(c('blue','lightgray','red')),
    strip = FALSE,at=at,
    colorkey = list(space='top')
  )
  lplots <- resizePanels(lplots,w=ngenes/nrow(dat))
  
  ann <- as.data.frame(
    sapply(
      more.genes,
      function(x) {
        sel <- sapply(x,function(y) which(dat$GeneID%in%y),simplify = F)
        res <- rep(NA,nrow(dat))
        mapply(
          function(sel,col) if(length(sel)>0) {res[sel] <<- more.gene.cols[col]},
          sel,names(sel)
        )
        return(res)
      # x <- as.data.frame(do.call(rbind,mapply(cbind,x,names(x))),stringsAsFactors=F)
      # x <- x[!duplicated(x[,1]),,drop=F]
      # row.names(x) <- x[,1]
      # x$col <- col[x[,2]]
      # x[row.names(dat$rna)[!row.names(dat$rna)%in%row.names(x)],] <- NA
      # return(x)
      },simplify = T
    ),
    row.names=dat$GeneID,
    stringsAsFactors=F
  )
  
  names.y <- 16
  dat$names.y <- names.y+1.5
  gene.names <- xyplot(
    names.y~GeneID|geneset,
    data = dat,
    layout=c(length(gene.sets),1),
    scales=list(
      x=list(
        rot=90,relation='free',cex=cex.axis,draw=F,axs='i'
      ),
      y=list(
        draw=T,cex=cex.axis,relation='free'
      )
    ),
    ylab=NULL,
    # yscale.component=yscale.components.default,
    between=list(x=2),                       #between creates more space
    # between the panels - you may need to adjust this value
    par.settings=list(layout.widths=list(right.padding=6)), #this creates
    yscale.component=function(...,right=F) yscale.components.default(...,right = F),
    panel = function(x, y,subscripts,...) {
      panel.xyplot(1:length(x),y,subscripts,...)
      panel.text(1:length(x),y,dat[subscripts,'gene'],srt=90,adj=c(1,.5))
      ann <- ann[subscripts,,drop=F]
      mapply(function(ann,yedge){
        # tmp <- ann$col
        # names(tmp) <- row.names(ann)
        mapply(
          panel.rect,
          1:length(x)-.5,
          names.y+yedge+.5,
          1:length(x)+.5,
          names.y+yedge+2.5,
          col=ann,
          border=NA
        )
      },ann,1:length(ann)*2)
    },
    prepanel=function(x, y, ...) {
      # x <- factor(x,x)
      return(list( 
        xlim = c(0,length(x)+1)#as.character(dat[levels(x),'gene']), 
        # xat = sort(unique(as.numeric(x)))
      ))
    },
    # prepanel=function(x, y, subscripts,...) {
    #   xsub <- unique(x[subscripts])
    #   return(list(
    #     xat = as.numeric(xsub)
    #   ))
    # },
    strip = FALSE,ylim=c(0,names.y+2.5+length(ann)*2),xlab=NULL
  )
  gene.names <- resizePanels(gene.names,w=ngenes/nrow(dat))
    
  # lplots <- resizePanels(lplots,w=ngenes/nrow(dat))
  height <- 3*(length(plots)+1.9)
  width <- sum(ngenes)/4+2
  # plots <- do.call(
  #   c,
  #   append(
  #     append(plots,list(gene.names),0), 
  #     list(
  #       lplots,
  #       merge.legends=F,
  #       layout=c(length(gene.sets),length(plots)+2)
  #     )
  #   )
  # )
  # plots <- resizePanels(plots,h=c(.90,rep(1,length(plots)-1)),w=ngenes/nrow(dat))
  
  # grid.newpage()
  # vp <- viewport(width = width/max(width,height),height = height/max(width,height))
  if(is.character(file)){
    dir.eps(file,path,width=height,height=width)
  }
  vp <- viewport(angle = 270,width = unit(width,'inches'),height = unit(height,'inches'))
  pushViewport(vp)
  
  i <- length(plots)+2.25
  print(lplots,position=c(0,(i-1.25)/i,1,1),more=T,newpage=F)
  mapply(function(x,y) print(x,position=c(0,(i-y-1.25)/i,1,(i-y-0.25)/i),more=T),plots,1:length(plots))
  print(gene.names,position=c(0,0,1,1/i))
  
  if(is.character(file)) dev.off()
  # popViewport()
  
  # plots <- append(
  #   append(plots,list(gene.names)), 
  #   list(
  #     lplots
  #   )
  # )
  # # do.call(grid.arrange,append(
  # #   plots,
  # #   list(layout_matrix=matrix(length(plots):1,length(plots)))
  # # ))
  # plots <- do.call(c,append(
  #   plots,
  #   list(
  #     merge.legends=F,
  #     layout=c(length(gene.sets),length(plots))
  #   )
  # ))
  # if(is.character(file)){
  #   dir.eps(file,path,width=width,height=height)
  #   plot(plots)
  #   dev.off()
  # }
  # 
  # pushViewport(viewport(unit(width,'inches'),unit(height,'inches')))
  # plot(plots)
  # rna <- lapply(rna,function(x) dbReadTable(con,x))
  # atac <- lapply(atac,function(x) dbReadTable(con,x))
  # gene.sets <- lapply(gene.sets,function(x) dbGetQuery(con,x))
  # peak.sets <- lapply(peak.sets,function(x) dbGetQuery(con,x))
}

diamondPlotPanel <- function(
  rna,atac,gene.peak,peak.sets,peak.cols,
  fn=panel.rna.atac,ylim='rna',p=0.05,cex.axis=1.2
){
  require(lattice)
  require(latticeExtra)
  
  atac <- sapply(
    Reduce(append,gene.peak),
    function(x) atac[x,'log2FoldChange',drop=F],
    simplify = F
  )
  
  dat <- as.data.frame(do.call(rbind,mapply(
    cbind,
    lapply(gene.peak,names),
    names(gene.peak),
    SIMPLIFY = F
  )),stringsAsFactors=F)
  names(dat) <- c("GeneID","geneset")
  dat <- dat[!duplicated(dat$GeneID),]
  dat$log2FoldChange <- rna[dat$GeneID,"log2FoldChange"]
  
  dat$sig <- rna[dat$GeneID,"padj"]<p&!is.na(rna[dat$GeneID,"padj"])
  # names(rnasig) <- dat$GeneID[rnasig]
  
  ngenes <-  unlist(table(dat$geneset))
  
  
  # if(ylim=='rna') {
  #   ylim <- range(dat$log2FoldChange)
  # }else {
  #   ylim <- range(unlist(atac))
  # }
  # ylim <- c(floor(ylim[1]),ceiling(ylim[2]))
  rna.ylim <- range(dat$log2FoldChange,na.rm = T)
  atac.ylim <- range(unlist(atac),na.rm = T)
  
  rna.atac.scale <- (atac.ylim[2]-atac.ylim[1])/(rna.ylim[2]-rna.ylim[1])
  atac <- lapply(atac,'/',rna.atac.scale)
  ylim <- range(c(rna.ylim,atac.ylim/rna.atac.scale))
  
  atac.tics <- yscale.components.default(
    range(c(atac.ylim,rna.ylim*rna.atac.scale))
  )
  #Custom y-scale component
  myyscale.component <- function(lim,...)
  {
    ans <- yscale.components.default(lim,...)
    # ans$right <- ans$left
    tmp <- yscale.components.default(lim*rna.atac.scale,...)#atac.tics$left
    ans$right <- tmp$left
    foo <- ans$right$labels$at
    # ans$right$labels$labels <- as.character(foo)
    ans$right$labels$at <- foo/rna.atac.scale
    ans$right$ticks$at <- foo/rna.atac.scale
    # scale factor difference between axes, adjust as necessary
    return(ans)
  }
  
  #The plot
  plots <- dotplot(
    log2FoldChange~GeneID|geneset,
    dat,
    layout=c(sum(ngenes>0),1),
    strip=FALSE,
    scales=list(x=list(
        rot=90,relation='free',cex=cex.axis,
        draw=F,axs="i"
    ),y=list(
      cex=cex.axis,draw=T,relation='free'
    )),
    # the panels
    yscale.component=myyscale.component,
    between=list(x=2),                       #between creates more space
    # between the panels - you may need to adjust this value
    par.settings=list(layout.widths=list(right.padding=6)), #this creates
    # more space on the right hand side of the plot
    ylab=NULL,
    xlab=NULL,
    prepanel=function(x, y, ...) {
      # x <- factor(x,x)
      return(list(
        xlim = c(0,length(x)+1),
        ylim=ylim
      ))
    },
    panel = function(x,y,groups,subscripts,...) {
      # print(subscripts)
      fn(
        x,y,groups,subscripts,...,atac=atac,
        sig=dat$sig,cex=2,
        peak.sets=peak.sets,peak.cols=peak.cols
      )
    }
  )
  plots <- resizePanels(plots,w=ngenes/nrow(dat))
  
  return(list(plots))
}
  
filtpts <- function(x,y,sel,...,cex=1.5){
  y <- y[sel,]
  x <- x[sel]
  # x <- x[sel]
  if(length(y)>0) return(panel.points(
    # rep(x,length(y))),
    x,y,...,cex=cex
  ))
}


panel.rna <- function(
  x,y,groups,subscripts,atac,
  cex.x=1.25,
  sig=T,cex=.8,
  peak.sets=NULL,peak.cols=c('red','blue'),...
){
  y[y=0] <- NA
  atac <- atac[subscripts]
  sig <- sig[subscripts]
  xval <- 1:length(x)
  panel.dotplot(xval,y,...,col=1,cex=cex.x,pch=NA,col.line = NA)
  panel.points(xval[!sig]-.25,y[!sig],pch=1,cex=cex.x,col=1)
  # if(!is.null(sig)){
    # sig <- sig[as.character(x)]
    panel.points(xval[sig]-.25,y[sig],pch=19,cex=cex.x,col=1)
  # }
  # draw lines at cutoffs
  # if(!is.null(peak.sets)) panel.abline(h=lfc,col=1)
  # draw line at 0
  panel.abline(h=0,col=1)
}

panel.atac <- function(
  x,y,groups,subscripts,atac,
  cex.x=1.25,
  sig=T,cex=.8,
  peak.sets=NULL,peak.cols=c('red','blue'),...
){
  y[y=0] <- NA
  atac <- atac[subscripts]
  sig <- sig[subscripts]
  xval <- 1:length(x)
  atac.x <- lapply(1:length(atac)+.25,function(z) jitter(rep(z,length(atac[[z]])),amount=.25))
  panel.dotplot(xval,y,...,col=1,cex=cex.x,pch=NA)
  mapply(filtpts,x=atac.x,y=atac,MoreArgs=list(
    sel=T,pch=18,col='gray',cex=cex
  ))
  peak.sel <- lapply(
    peak.sets,
    function(i) lapply(atac, function(j) row.names(j)%in%i)
  )
  mapply(function(peak.sel,peak.cols) mapply(
    filtpts,x=atac.x,y=atac,sel=peak.sel,
    MoreArgs = list(pch=18,cex=cex,col=peak.cols)
  ),peak.sel,peak.cols)
  # panel.points(xval[!sig]-.25,y[!sig],pch=1,cex=cex.x,col=1)
  # panel.points(xval[sig]-.25,y[sig],pch=19,cex=cex.x,col=1)
  # draw lines at cutoffs
  # if(!is.null(peak.sets)) panel.abline(h=lfc,col=1)
  # draw line at 0
  panel.abline(h=0,col=1,lty = 2)
}

panel.rna.atac <- function(
  x,y,groups,subscripts,atac,
  cex.x=1.25,
  sig=T,cex=.8,
  peak.sets=NULL,peak.cols=c('red','blue'),...
){
  y[y=0] <- NA
  atac <- atac[subscripts]
  sig <- sig[subscripts]
  xval <- 1:length(x)
  atac.x <- lapply(1:length(atac)+.25,function(z) jitter(rep(z,nrow(atac[[z]])),amount=.25))
  panel.dotplot(xval,y,...,col=1,cex=cex.x,pch=NA)
  mapply(filtpts,x=atac.x,y=atac,MoreArgs=list(
    sel=T,pch=18,col='gray',cex=cex
  ))
  peak.sel <- lapply(
    peak.sets,
    function(i) lapply(atac, function(j) row.names(j)%in%i)
  )
  mapply(function(peak.sel,peak.cols) mapply(
    filtpts,x=atac.x,y=atac,sel=peak.sel,
    MoreArgs = list(pch=18,cex=cex,col=peak.cols)
  ),peak.sel,peak.cols)
  panel.points(xval[!sig]-.25,y[!sig],pch=1,cex=cex.x,col=1)
  panel.points(xval[sig]-.25,y[sig],pch=19,cex=cex.x,col=1)
  # draw lines at cutoffs
  # if(!is.null(peak.sets)) panel.abline(h=lfc,col=1)
  # draw line at 0
  panel.abline(h=0,col=1,lty = 2)
}

