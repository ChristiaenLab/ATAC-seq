source('DESeqFns.R')

mergeRnaAtac <- function(
  rna,atac,genes,ann,gene.names,
  p=c(.05,.05),lfc=c(0,0),rnasigfn=is.sig,which.sig=1:length(rna),
  mask=F,df.fn=NULL,rmdup=T,fromLast=F,sig.sub=T,anysigfn=any,notsig=F
){
    lapply(1:length(genes),function(x) {
      if(!is.data.frame(genes[[x]])) genes[[x]] <<- data.frame(
        Type=rep('Uncertain',length(genes[[x]])),
        row.names = genes[[x]]
      ) else genes[[x]] <<- genes[[x]][,'Type',drop=F]
      if(!"Type"%in%names(genes[[x]])) genes[[x]]$Type <<- "Uncertain"
      genes[[x]]$split <<- x
    })
    dfs <- genes
    genes <- sapply(genes,row.names,simplify = F)
  if(rmdup){
    dup <- unlist(genes)[duplicated(unlist(genes))]
    genes <- lapply(genes,function(x) x[!x%in%dup])
  }
  dfs <- mapply(function(x,y) x[y,,drop=F],dfs,genes)
  split <- factor(
    unlist(dfs['split',]),
    levels = names(genes)
  )
  rna <- lapply(rna,'[',as.character(unlist(genes)),)
  sig <- mapply(rnasigfn,rna, p=p[[1]],lfc=lfc[[1]])
  sig[is.na(sig)] <- F
  row.names(sig) <- as.character(unlist(genes))
  rnasel <- data.frame(
    split,
    sig=T,
    Type=unlist(dfs['Type',]),
    genes = unlist(genes),stringsAsFactors = F
  )
  if(!is.null(which.sig)){
    rnasel$sig <- apply(sig[,which.sig,drop=F],1,anysigfn)
  }
  rnasel <- subset(rnasel,!duplicated(genes,fromLast=fromLast))
  row.names(rnasel) <- rnasel$genes
  rna <- lapply(rna,'[',row.names(rnasel),)
  sig <- sig[row.names(rnasel),,drop=F]
  # sig <- sapply(rna, function(x) x$padj<p[1]&abs(x$log2FoldChange)>lfc[1])
  # row.names(sig) <- row.names(rnasel)
  if(mask) {
    rna <- mapply(function(x,sel){
      x[sel,'log2FoldChange'] <- NA
      return(x[,'log2FoldChange'])
    },x=rna,sel=as.list(sig)) 
    rnasel$na <- T
  }else {
    rna <- sapply(rna,'[',,'log2FoldChange')
    rna[is.na(rna)] <- 0
  }
  row.names(rna) <- row.names(rnasel)
  
  if(sig.sub) rnasel <- subset(rnasel,sig)
  rna <- rna[row.names(rnasel),,drop=F]
  sig <- sig[row.names(rnasel),,drop=F]
  
  atacsel <- row.names(atac[[1]])[any.sig(atac,lfc=lfc[[2]],p=p[[2]])]
  atac <- sapply(atac,'[',atacsel,'log2FoldChange')
  row.names(atac) <- atacsel
  
  peaks <- ann$geneToPeak[row.names(rnasel)]
  peaks <- sapply(peaks,function(x) Filter(function(y) y%in%row.names(atac),x))
  peaks <- Filter(function(x) length(x)>0,peaks)
  peaks <- lapply(peaks,function(y) atac[y,,drop=F])
  
  rna <- rna[names(peaks),,drop=F]
  rnasel <- rnasel[names(peaks),,drop=F]
  sig <- sig[names(peaks),,drop=F]
  rnasel$gene <- gene.names[row.names(rna),]
  rnasel$gene <- sub('KH2013:','',rnasel$gene)
  rnasel$gene <- sub('.*_','',rnasel$gene)
  rnasel$gene <- sub('(.{,12}).*','\\1',rnasel$gene)
  # rnasel$gene <- sub('/.*','',rnasel$gene)
  return(list(rna=rna,atac=atac,peaks=peaks,rnasel=rnasel,sig=sig))
}

filtpts <- function(fn,x,y,...,cex=1.5){
  sel <- fn(y)
  y <- y[sel]
  x <- x[sel]
  if(length(y)>0) return(panel.points(
    # rep(x,length(y))),
    x,y,...,cex=cex
  ))
}
panel.rnaAtac.dot <- function(
  x,y,atac,labels,lfc,ylim,
  cex.x=1.25,
  ann=NULL,ma=NULL,sig=NULL,cex=.8,
  ref.sig=NULL,peakcols=c('red','blue'),...
){
  x <- factor(x,x)
  atac <- atac[names(y)]
  # panel.abline(h=seq(ylim[1],ylim[2],.5),col='gray')
  if(!is.null(ref.sig)){
    filters <- lapply(
      ref.sig,function(x) return(
        function(w) relist(row.names(w)%in%x,w)
      )
    )
  }else{
    filters <- list(
      function(w) w>lfc[2],
      function(w) w< lfc[1]
    )
  }
  filt.gray <- function(w) !Reduce(
    function(i,j) i|j,
    lapply(filters,function(fn) fn(w))
  )
  panel.dotplot(x,y,...,col=1,cex=cex.x,pch=NA)
  atac.x <- lapply(1:length(atac),function(z) jitter(rep(z,length(atac[[z]])),amount=.25))
  mapply(filtpts,x=atac.x,y=atac,MoreArgs=list(
    fn=filt.gray,pch=18,col='gray',cex=cex
  ))
  # lapply(1:length(atac),function(z){
  #   filtpts(filt.gray,atac.x[[z]],atac[[z]],pch=18,col='gray',cex=cex)
  # })
  mapply(function(x,y) mapply(
    filtpts,fn=filters,col=peakcols,
    MoreArgs = list(x=x,y=y,pch=18,cex=cex)
    # filtpts(filt.red,z,atac[[z]],pch=18,col='red',cex=cex)
    # filtpts(filt.blue,z,atac[[z]],pch=18,col='blue',cex=cex)
  ),x=atac.x,y=atac)
  panel.points(x[y!=0],y[y!=0],pch=1,cex=cex.x,col=1)
  if(!is.null(sig)){
    sig <- sig[as.character(x)]
    panel.points(x[sig],y[sig],pch=19,cex=cex.x,col=1)
  }
  # draw lines at cutoffs
  # if(!is.null(ref.sig)) panel.abline(h=lfc,col=1)
  # draw line at 0
  panel.abline(h=0,col=1)
  if(!is.null(ann)) {
    tmp <- ann$col
    names(tmp) <- row.names(ann)
    mapply(
      panel.rect,
      1:length(x)-.5,
      ylim[1]-.5,
      1:length(x)+.5,
      ylim[1],
      col=ann[as.character(x),'col'],
      border=NA
    )
  }
}

rnaAtacPlotPanel <- function(
  rna,atac,labels,file=F,path='.',lfc=c(-1,1),sig=NULL,
  ...,plot='all',order=T,ann=NULL,ma=NULL,ref.sig=NULL,
  peakcols=c('blue','red'),cex.axis=.6
){
  require(lattice)
  require(latticeExtra)
  names(rna) <- row.names(labels)
  if(!is.null(sig)) names(sig) <- row.names(labels)
  if(plot!='all') {
    sel <- sapply(atac,function(x) any(x<lfc[1]|x>lfc[2]))
    if(!is.na(ref.sig)) sel <- sel&sapply(
      atac,function(x) any(row.names(x)%in%unlist(ref.sig))
    )
    if(plot=='notDA') sel <- !sel
    rna <- rna[names(atac)[sel]]
    rna <- na.omit(rna)
  }
  if(order) rna <- rna[order(rna)]
  atac <- atac[names(rna)]
  labels <- labels[names(rna),]
  labels <- do.call(
    data.frame,
    append(labels,list(stringsAsFactors=T,row.names=row.names(labels)))
  )
  ylim <- range(c(rna,unlist(atac)))
  ylim <- c(floor(ylim[1]),ceiling(ylim[2]))
  ngenes <- as.vector(table(labels[,c("Type",'split')]))
  
  key <- list(
    columns=4,
    space='bottom',
    points=list(
      pch=c(1,19,18,18,18),
      col=c('black','black','gray','blue','red')
    ),
    text=list(lab=c(
      "not significant",
      "DE gene (FDR < 0.05)",
      "DA peak (FDR < 0.05)",
      paste0("log2FC < ",as.character(lfc[1])),
      paste0("log2FC > ",as.character(lfc[2]))
    ))
  )
  legend <- list(top=list(fun=draw.key,args=list(key=key)))
  
  if(!is.null(ann)) {
    legend$bottom <- list(
      fun=draw.key,
      args=list(key=list(
        columns=length(unique(ann$col)),
        rectangles=list(col=unique(ann$col)),
        text=list(lab=unique(ann[,2]))
      ))
    )
    key$rectangles <- list(col=unique(ann$col))
    key$text$lab <- c(key$text$lab,unique(ann[,2]))
  }
  plots <- dotplot(
    rna~factor(names(rna),names(rna))|labels[,'Type']+labels[,'split'],
    layout=c(length(unique(labels$split))*length(levels(labels[,"Type"])),1),
    strip=FALSE,
    scales=list(x=list(
        rot=90,relation='free',cex=cex.axis,draw=F
    ),y=list(cex=cex.axis,draw=T)),ylab='log2FC',
    prepanel=function(x, y, ...) {
      x <- factor(x,x)
      return(list( 
        xlim = as.character(labels[levels(x),'gene']), 
        xat = sort(unique(as.numeric(x))),ylim=ylim
      ))
    },
    panel = function(...) panel.rnaAtac.dot(
      ...,atac=atac,labels = labels,ylim = ylim,lfc=lfc,ann=NULL,
      ma=ma,sig=sig,cex=2,ref.sig=ref.sig,peakcols=peakcols
    ),
    legend=legend
  )
  plots <- resizePanels(plots,w=ngenes/length(rna))
  
  if(is.character(file)){
    dir.eps(file,path,width=16,height=8)
    plot(plots)
    dev.off()
    dir.csv(data.frame(
      rnaLFC=rna,
      GeneName=labels$gene,group=labels$split,type=labels$Type,
      peaks=sapply(atac,function(x) paste(row.names(x),collapse=';')),
      atacLFC=sapply(atac,function(x) paste(as.character(x),collapse=';'))
    ),file,path)
  }
  return(list(plots))
}
  
rnaAtacPanelHmap <- function(
  rna,atac,genes,ann,gene.names,file,ma,
  p=c(.05,.05),lfc=c(0,0),rnasigfn=is.sig,which.sig=1:length(rna),
  path='.',cex.axis=1,col.lfc=0,levels=NULL,col=NULL,which.da.lib=NULL,which.da.gene=any,
  width=NULL,height=NULL,sig.sub=T,anysigfn=any,ma.rm=F,dend=T,
  ref.sig=NULL,peakcols=c('blue','red'),...
){
  col.lfc <- lapply(lfc[[2]],function(i) if(length(i)==1) {c(-i,i)} else i)
  # if(!is.null(dafilt)){
  #   sel <- mapply(is.sig,atac,lfc[[2]],p=p[[2]],SIMPLIFY = F,USE.NAMES = F)
  #   # sel <- sapply(
  #   #   1:length(col.lfc),
  #   #   function(i) sapply(
  #   #     dat$peaks,
  #   #     function(x) any(x[,i]<col.lfc[[i]][1]|x[,i]>col.lfc[[i]][2])
  #   #   )
  #   # )
  #   sel <- do.call(dafilt,sel)
  #   atacsig <- row.names(atac[[1]])[sel]
  #   #uncomment to remove nonsignificant peaks 
  #   # atac <- sapply(atac,'[',sel,T)
  #   
  #   # dat$peaks <- dat$peaks[sel]
  #   # dat$rna <- dat$rna[sel,]
  #   # dat$rnasel <- dat$rnasel[sel,]
  #   # dat$sig <- dat$sig[sel,]
  # }
  dat <- mergeRnaAtac(
    rna,atac,genes,ann,gene.names,p=p,lfc=lfc,rnasigfn,which.sig=which.sig,mask = F,
    sig.sub = sig.sub,anysigfn = anysigfn
  )
  if(!is.null(which.da.lib)){
    sel <- mapply(is.sig,atac,lfc[[2]],p=p[[2]],SIMPLIFY = F,USE.NAMES = F)
    sel <- do.call(which.da.lib,sel)
    atacsig <- row.names(atac[[1]])[sel]
    # sel <- sapply(
    #   1:length(col.lfc),
    #   function(i) sapply(
    #     dat$peaks,
    #     function(x) any(x[,i]<col.lfc[[i]][1]|x[,i]>col.lfc[[i]][2])
    #   ),simplify = F,USE.NAMES = F
    # )
    # sel <- do.call(dafilt,sel)
    sel <- sapply(dat$peaks,function(x) which.da.gene(row.names(x)%in%atacsig))
    dat$peaks <- dat$peaks[sel]
    dat$rna <- dat$rna[sel,]
    dat$rnasel <- dat$rnasel[sel,]
    dat$sig <- dat$sig[sel,]
  }
  if(!is.null(levels)&!is.null(col)){
    if(!is.list(levels)) levels <- rep(levels,length(rna))
    levels <- lapply(levels,function(x) {
      x <- as.data.frame(do.call(rbind,mapply(cbind,x,names(x))),stringsAsFactors=F)
      x <- x[!duplicated(x[,1]),,drop=F]
      row.names(x) <- x[,1]
      x$col <- col[x[,2]]
      x[row.names(dat$rna)[!row.names(dat$rna)%in%row.names(x)],] <- NA
      return(x)
    })
  }
  if(ma.rm){
    sel <- intersect(row.names(dat$rna),row.names(ma))
    ma <- ma[sel,]
    dat$peaks <- dat$peaks[sel]
    dat$rna <- dat$rna[sel,]
    dat$rnasel <- dat$rnasel[sel,]
    dat$sig <- dat$sig[sel,]
  }else{
    sel <- row.names(dat$rna)[!row.names(dat$rna)%in%row.names(ma)]
    ma[sel,] <- 0
    ma <- ma[row.names(dat$rna),]
  }
  if(dend){
    dend <- as.dendrogram(hclust(dist(as.matrix(ma))))
    sel <- row.names(ma)[order.dendrogram(dend)]
  }else sel <- row.names(dat$rna)[order(dat$rna[,min(which.sig)])]
  sel <- do.call(c,split(sel,paste0(
    dat$rnasel[sel,'split'],
    as.character(as.numeric(as.factor(dat$rnasel[sel,'Type'])))
  )))
  ma <- ma[sel,,drop=F]
  rna <- dat$rna[sel,,drop=F]
  peaks <- dat$peaks[sel,drop=F]
  rnasel <- dat$rnasel[sel,,drop=F]
  sig <- dat$sig[sel,,drop=F]
  ngenes <- as.vector(table(rnasel[,c("Type",'split')]))
  ngroups <- length(unique(rnasel$split))*length(unique(rnasel[,"Type"]))
  plots <- mapply(
    rnaAtacPlotPanel,as.data.frame(rna),
    lapply(names(atac),function(x) lapply(peaks,'[',,x,drop=F)),
    path=paste0(path,'/',names(rna)),
    lfc=col.lfc,
    sig=as.data.frame(sig),#ann=levels,
    MoreArgs = list(
      labels=rnasel,order=F,ref.sig=ref.sig,peakcols=peakcols,cex.axis=cex.axis,...
    )
  )
  zmax <- quantile(abs(unlist(ma)),.99)
  ma[ma>zmax] <- zmax
  ma[ma< -zmax] <- -zmax
  at <- seq(-zmax,zmax,length.out = 16)
  lplots <- levelplot(
    z~x*y|type*split,
    data = data.frame(
      expand.grid(
        x=factor(row.names(ma),row.names(ma)),
        y=names(ma)
      ),
      z=unlist(ma),
      split=rnasel[,'split'],
      type=rnasel[,'Type']
    ),
    layout=c(ngroups,1),
    scales=list(
      x=list(
        rot=90,relation='free',cex=0,draw=F
      ),
      y=list(cex=cex.axis,draw=T)),
    prepanel=function(x, y, subscripts,...) {
      xsub <- unique(x[subscripts])
      return(list( 
        xat = as.numeric(xsub)
      ))
    },
    col.regions = colorRampPalette(c('blue','lightgray','red')),
    strip = FALSE,at=at
  )
  lplots <- resizePanels(lplots,w=ngenes/length(rna))
  
  names.y <- 16
  gene.names <- xyplot(
    rep(names.y+1.5,nrow(rnasel))~genes|split,
    data = rnasel,
    layout=c(ngroups,3),
    scales=list(x=list(
        rot=90,relation='free',cex=cex.axis,draw=F
    ),y=list(draw=T)),
    panel = function(x, y,...,ann=levels) {
      panel.xyplot(x,y,...)
      panel.text(1:length(x),y,rnasel[as.character(x),'gene'],srt=90,adj=c(1,.5))
      if(!is.null(ann)) mapply(function(ann,yedge){
        tmp <- ann$col
        names(tmp) <- row.names(ann)
        mapply(
          panel.rect,
          1:length(x)-.5,
          names.y+yedge+.5,
          1:length(x)+.5,
          names.y+yedge+2.5,
          col=ann[as.character(x),'col'],
          border=NA
        )
      },ann,1:length(ann)*2)
    },
    prepanel=function(x, y, ...) {
      x <- factor(x,x)
      return(list( 
        xlim = as.character(rnasel[levels(x),'gene']), 
        xat = sort(unique(as.numeric(x)))
      ))
    },
    # prepanel=function(x, y, subscripts,...) {
    #   xsub <- unique(x[subscripts])
    #   return(list(
    #     xat = as.numeric(xsub)
    #   ))
    # },
    strip = FALSE,ylim=c(0,names.y+2.5+length(levels)*2),xlab='',ylab=''
  )
  gene.names <- resizePanels(gene.names,w=ngenes/length(rna))
    
  lplots <- resizePanels(lplots,w=ngenes/length(rna))
  if(is.null(height)) height <- 2*(length(plots)+1.9)
  if(is.null(width)) width <- nrow(rna)/4
  plots <- Reduce(
    function(x,y) c(x,y,merge.legends=F,layout=c(ngroups,length(plots)+2)),
    append(append(plots,list(gene.names),0), list( lplots )
  ))
  plots <- resizePanels(plots,h=c(.90,rep(1,length(plots)-1)),w=ngenes/length(rna))
  dir.pdf(file,path,width=width,height=height)
  plot(plots)
  dev.off()
}

