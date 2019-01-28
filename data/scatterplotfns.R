lscatterhist <- function(
  x,y,xlab='',ylab='',legend=NULL,lpos='bottomright',
  pch=1:length(x),col=1:length(x),pnames=NULL,
  groups=NULL,histlegend=NULL,#main=F,
  ylim=NA,xlim=NA,cex=1,which.hists=1:length(x),fit=NULL,
  pcex=cex,cutoff=NULL,segments=F,...){
  x <- lapply(x,unlist)
  y <- lapply(y,unlist)
  sapply(1:length(x),function(i) x[[i]][!is.finite(x[[i]])] <<- NA)
  sapply(1:length(y),function(i) y[[i]][!is.finite(y[[i]])] <<- NA)
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  if(any(which.hists)) layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  if(is.na(xlim)) xlim <- range(unlist(x),na.rm = T)
  if(is.na(ylim)) ylim <- range(unlist(y),na.rm = T)
  sapply(1:length(x),function(i) {
    sel <- x[[i]]>=xlim[1]&x[[i]]<=xlim[2]&y[[i]]>=ylim[1]&y[[i]]<=ylim[2]
    x[[i]] <<- x[[i]][sel]
    y[[i]] <<- y[[i]][sel]
    # if(length(pnames[[i]])>0) pnames[[i]] <<- pnames[[i]][sel]
  })
    
  if(length(groups)==length(x)){
    xhist = lapply(na.omit(unique(groups)),function(i) hist(
      unlist(x[groups==i]), plot=FALSE,
      breaks = seq(floor(xlim[1]),ceiling(xlim[2]),length.out = 30)))
    yhist = lapply(na.omit(unique(groups)),function(i) hist(
      unlist(y[groups==i]), plot=FALSE,
      breaks = seq(floor(ylim[1]),ceiling(ylim[2]),length.out = 30)))
    # top = max(Reduce(function(i,j) c(i, xhist[[j]]$counts,yhist[[j]]$counts),1:length(xhist)))
    histcol <- sapply(na.omit(unique(groups)),function(i) col[groups==i][1])
  }
  else{
    xhist = lapply(x,function(i) hist(
      i, plot=FALSE,
      breaks = seq(floor(xlim[1]),ceiling(xlim[2]),length.out = 30)))
    yhist = lapply(y,function(i) hist(
      i, plot=FALSE,breaks = seq(floor(ylim[1]),ceiling(ylim[2]),length.out = 30)))
      # function(i,j) c(i, xhist[[j]]$counts,yhist[[j]]$counts),1:length(xhist)))
    histcol <- col
    # histcol=sapply(col,adjustcolor,.5)
  }
  xhist <- do.call(rbind,lapply(xhist,function(x) x$counts/sum(x$counts,na.rm = T)))
  yhist <- do.call(rbind,lapply(yhist,function(x) x$counts/sum(x$counts,na.rm = T)))
  topx <- apply(xhist,2,sum)
  topy <- apply(yhist,2,sum)
  topx = max(unlist(topx),na.rm = T)
  topy = max(unlist(topy),na.rm = T)
  par(mar=c(4,4,2,1))
  # if(main) main <- paste(
  #   'r =',as.character(round(cor.test(x[[1]],y[[1]])$estimate,2)),
  #   '\trho =',as.character(round(cor.test(x[[1]],y[[1]],method = 'spearman')$estimate,2))
  #   )else main <- ''
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,
       ...)
  # if(is.character(legend)) legend(lpos,pch = pch[which.hists],col=col[which.hists],legend[which.hists])
  mapply(points,x,y,col=col,pch=pch,cex=pcex)
  if(sum(fit)>0) {
    sapply(which(fit), function(i) abline(
      lm(y[[i]]~x[[i]]),col=col[i]
    ))
    # abline(
      # lm(y~0+x,data=data.frame(x=unlist(x[fit]),y=unlist(y[fit]))
    # ))
  }
  if(length(pnames)>0){
    # sapply(
    # 1:length(pnames),
    # function(i) if(length(pnames[[i]])>0) {
    text(
      pnames$x,pnames$y,pnames$labs,cex = cex/2,adj=c(.5,0)
    )
    # )}
  }
  if(length(cutoff)==2) {
    if(segments) segments(
      c(rep(xlim[1],2),-cutoff[1],cutoff[1]),
      c(-cutoff[2],cutoff[2],ylim[1],ylim[1]),
      c(rep(xlim[2],2),-cutoff[1],cutoff[1]),
      c(-cutoff[2],cutoff[2],ylim[2],ylim[2])
    ) else {
      # rect(xlim[1],-cutoff[2],xlim[2],cutoff[2],col = rgb(1,1,1,.5),border = NA)
      # rect(-cutoff[2],ylim[1],-cutoff[2],ylim[2],col = rgb(1,1,1,.5),border = NA)
    }
  }
  par(mar=c(0,3,1,1))
  # add <- F
  # sapply(which.hists,function(i) {
  #   barplot(
  #     xhist[[i]], 
  #     axes=FALSE, ylim=c(0, topx), space=0,col=histcol[i],add=add)
  #   add <<- T})
  if(any(which.hists)){
    barplot(
      xhist[which.hists,],
      beside = F,
      col = histcol[which.hists],
      border = NA,ylim=c(0, topx),space = 0,ylab="Density")
    legend('topright',fill=col[which.hists],legend[which.hists],cex=.8)
    par(mar=c(3,0,1,1))
    # add <- F
    # sapply(which.hists,function(i) {
    #   barplot(
    #     yhist[[i]], 
    #     axes=FALSE, xlim=c(0, topy), space=0, horiz=TRUE,
    #     col = histcol[i],add=add)
    #   add <<- T})
    barplot(
      yhist[which.hists,],
      beside = F,
      col = histcol[which.hists],
      border = NA,xlim=c(0, topy),space = 0,horiz = T,xlab = 'Density')
    par(oma=c(3,3,0,0))
    # mtext(xlab, side=1, line=2.5, outer=TRUE, adj=0, 
    #       at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
    # mtext(ylab, side=2, line=2, outer=TRUE, adj=0, 
    #       at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
  }
}

lbarplot <- function(x,legend=NULL,col=1:length(x),lpos='topleft',xlim=NULL,...){
  sapply(1:length(x),function(i) x[[i]][!is.finite(x[[i]])] <<- NA)
  if(is.null(xlim)) xlim=range(unlist(x),na.rm = T)
  sapply(1:length(x),function(i) {
    sel <- x[[i]]>=xlim[1]&x[[i]]<=xlim[2]
    x[[i]] <<- x[[i]][sel]
  })
  breaks <- seq(floor(xlim[1]),ceiling(xlim[2]),length.out = 30)
  xhist = lapply(x,function(i) hist(
    i, plot=FALSE,breaks = breaks
    ))
  xhist <- do.call(rbind,lapply(xhist,function(i) i$counts/sum(i$counts,na.rm = T)))
  # histcol=sapply(col,adjustcolor,.5)
  barplot(
    xhist,
    beside = F,
    col = col,
    border = NA,space = 0,ylab='Density')
  axis(1,c(0,5,10,15,20,25,30),labels = c(0,round(breaks)[c(5,10,15,20,25,30)]))
  
  # add <- F
  # sapply(1:length(xhist),function(i) {
  #   barplot(
  #     xhist[[i]]$counts/sum(xhist[[i]]$counts,na.rm = T), 
  #     col=histcol[i],add=add,...)
  #   add <<- T})
  if(is.character(legend)) legend(lpos,fill=col,legend=legend)
}

cor.title <- function(x,y) paste(
  'r =',as.character(round(cor(x,y),4)),
  'rho =',as.character(round(cor(x,y,method = 'spearman'),4))
)


diamondRnaAtac <- function(
  rna,atac,gene.names,filename,
  path='output',npeak=F,p=.05,lfc=1,
  ylab='log2(RNAseq FC)',...){
  rna <- rna[rna$padj<p&abs(rna$log2FoldChange)>lfc,]
  # sel <- grep('KH2013:',labs)
  # labs <- labs[-sel]
  # rna <- rna[-sel,]
  z <- geneToPeak(
    row.names(rna),
    atac)
  z <- lapply(z,function(i) atac[i,])
  if(npeak){
    z <- lapply(z,function(i) i[order(i$pvalue),'log2FoldChange'])
    z <- lapply(z,function(i) i[1:min(npeak,length(i))])
  }else{
    z <- lapply(z,function(i) i[i$padj<p,'log2FoldChange'])
    z <- lapply(z,function(i) i[order(i)])
  }
  z <- z[sapply(z,length)>0]
  # sapply(1:length(z),function(i) if(length(z[[i]])>npeak) z[[i]] <- z[[i]][1:npeak])
  rna <- rna[names(z),]
  labs <- gene.names[row.names(rna),]
  labs <- sub('KH2013:','',labs)
  labs <- sub('.*_','',labs)
  
  out <- data.frame(GeneID=row.names(rna),GeneName=labs,logFC=rna$log2FoldChange,npeaks=sapply(z,length))
  out <- out[order(out$logFC),]
  dir.tab(out,filename,path,row.names=F)
  dir.png(filename,path)
  diamondplot(rna$log2FoldChange,z,labs,ylab=ylab,...)
  dev.off()
}


diamondplot <- function(
  y,z,labs='',zlim=500,#filename, path='~/ciona',
  col=rev(brewer.pal(n=11, name="RdBu")),
  ysep=.05,pch=18,cex=80,...){
  require(RColorBrewer)
  require(fields)
  col <- colorRamp(col)
  getcol <- function(z){
    z <- z-min(z)
    z <- z/max(z)
    z <- col(z)
    z <- apply(
      z,1, 
      function(i) do.call(
        function(...) rgb(...,maxColorValue = 255),
        as.list(i)))
    return(z)
  }
  # plot(xlim=range(x),ylim = range(y),...)
  # sort data
  sel <- sapply(z,length)<zlim
  y <- y[sel]
  z <- z[sel]
  labs <- labs[sel]
  x <- seq(0,by=1,length.out = length(y))
  sel <- order(y)
  labs <- labs[sel]
  y <- y[sel]
  z <- z[sel]
  z <- lapply(z,sort)
  xlabs <- x
  ylabs <- y
  
  ylim <- range(unlist(y))
  ysep <- abs(ylim[2]-ylim[1])*ysep
  y <- lapply(1:length(y),function(i) seq(y[i],by=ysep,length.out=length(z[[i]])))
  x <- lapply(1:length(x),function(i) rep(x[i],length(z[[i]])))
  
  x <- unlist(x)
  y <- unlist(y)
  z <- unlist(z)
  color <- getcol(z)
  # dir.eps(filename,path)
  # zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  # layout(c(0,1), widths=c(4/5,1/5))
  ylim <- range(y)
  
  ylim <- c(
    min(ylabs-(nchar(labs)*(1/length(xlabs))*(ylim[2]-ylim[1]))),
    max(y)
  )
  plot(x,y,col=color,pch=pch,xaxt='n',xlab='',ylim=ylim,cex=cex/length(xlabs),...)
  text(
    xlabs,#+150/length(xlabs),
    ylabs,labs,cex = (cex/2)/length(xlabs),srt=90,
    # pos = 2,offset = 1
    adj=c(1,.5)
  )
  # colorbar.plot(
  #   max(x),min(y),round(seq(min(z),max(z),length.out = 11),2),
  #   col = rev(brewer.pal(11,'RdBu')),adj.x = 1,adj.y = 1,horizontal = F)
  # dev.off()
}

atacBoxplot <- function(atac,peaklist,fdr=0.05){
  atac <- subset(atac,padj<fdr)
  res <- lapply(peaklist,function(x) atac[x,'log2FoldChange'])
  boxplot(res,las=1)
}
latacBoxplot <- function(ataclist,peaklist,col=0,fdr=0.05,...){
  ataclist <- lapply(ataclist,function(x) x[x$padj<fdr,])
  res <- unlist(lapply(
    ataclist,
    function(atac) lapply(peaklist,function(x) atac[x,'log2FoldChange'])
  ),F)
  if(length(col)==length(peaklist)) col <- rep(col,length(ataclist))
  boxplot(res,col=col,las=2,...)
}

scatterGS <- function(
  x,y,genes=NULL,filename,path=gsub('\\W','_',paste(lab,collapse = '_')),
  gene.names=NULL,
  lab=names(genes),
  col=c('firebrick','chartreuse4'),plot.sig=T,plot.all=F,
  pch=rep(19,length(genes)),
  fit=F,total.fit=F,lfc=c(1,1),cex=rep(1.5,length(genes)+1),
  fdr=c(0.05,0.05),...
){
  dat <- merge(x,y,0)
  row.names(dat) <- dat$Row.names
  sig <- subset(dat,padj.x<fdr[1]|padj.y<fdr[2])
  de <- row.names(sig)[abs(sig$log2FoldChange.x)>lfc[1]&abs(sig$log2FoldChange.y)>lfc[2]]
  notDE <- row.names(sig)[!row.names(sig)%in%de]
  sel <- do.call(append,lapply(
    list(de,notDE),
    function(z) lapply( genes,intersect,z )
  ))
  cex <- c(cex,cex[-1])
  pch <- c(pch,pch[-1])
  col <- c(col,sapply(col,adjustcolor,offset=c(0,.6,.6,.6)))
  xyarg <- lapply(
      sel, function(w) sig[which(sig$Row.names%in%w),,drop=F]
    )
  do.call(lscatter,list(
    sig,sel,filename,path,gene.names,lab,
    col,plot.all=plot.sig,pch=pch,cex=cex,...
  ))
}

scatterRnaAtac <- function(
  rna,atac,ann,genes=NULL,filename,path=gsub('\\W','_',paste(lab,collapse = '_')),
  gene.names=NULL, geneCols=c('GeneID','GeneID10kUp','GeneID10kDown'),
  lfc=c(1,1),fdr=c(0.05,0.05),DE=T,...
){
  if(DE) genes <- append(genes,list(DA=row.names(sig.sub(rna,lfc[1],fdr[1]))),0)
  dat <- mergeRnaAtac(
    list(rna),list(atac),genes,ann,gene.names,
    rmdup=F,fromLast=T,lfc=list(lfc[1],0),p=c(fdr[1],1)
  )
  trpkdat <- do.call(rbind,mapply(
    cbind,
    row.names(dat$rna),
    sapply(dat$peaks,row.names),
    dat$rna,dat$peaks
  ))
  trpkdat <- as.data.frame(trpkdat,stringsAsFactors = F,row.names = F)
  names(trpkdat) <- c("Row.names",'peak',"log2FoldChange.x","log2FoldChange.y")
  sapply(3:4, function(x) trpkdat[,x] <<- as.numeric(trpkdat[,x]))
  sig <- trpkdat[!apply(trpkdat,1,anyNA),]
  sig <- sig[atac[sig[,2],'padj']<fdr[2],]
  vsig <- sig[abs(sig[,3])>lfc[1]&abs(sig[,4])>lfc[2],]
  de <- sig$Row.names[abs(sig$log2FoldChange.x)>lfc[1]&abs(sig$log2FoldChange.y)>lfc[2]]
  sel <- lapply(genes,intersect,de)
  xyarg <- lapply(
      sel, function(w) vsig[which(vsig$Row.names%in%w),,drop=F]
    )
  do.call(lscatter,list(
    trpkdat,xyarg,sel,filename,path,gene.names,...
  ))
}

lscatter <- function(
  dat,xyarg,sel=NULL,filename,path=gsub('\\W','_',paste(lab,collapse = '_')),
  gene.names=NULL,
  lab=names(sel),
  col=c('firebrick','chartreuse4'),plot.all=T,
  pch=rep(19,length(sel)+1),
  cex=rep(2,length(sel)+1),legend=c('FDR < 0.05',names(sel)),
  ...
){
  sel <- sapply(xyarg,nrow)>0
  xyarg <- xyarg[sel]
  col <- col[sel]
  pcex <- cex[c(F,sel)]
  ppch <- pch[c(F,sel)]
  legend <- legend[c(T,sel)]
  pnames <- NULL
  if(!is.null(gene.names)) {
    labs <- do.call(rbind,lapply(
      xyarg,'[',,c('Row.names','log2FoldChange.x','log2FoldChange.y')
    ))
    labs <- labs[order(abs(labs$log2FoldChange.y),decreasing = T),]
    labs <- labs[!duplicated(labs$Row.names),]
    pnames <- data.frame(
      labels=sub('KH2013:','',sub(
        'KH2013:.*_','',sub('/.*','',gene.names[labs[,1],])
      )),
      x=as.numeric(labs[,2]),
      y=as.numeric(labs[,3]),
      stringsAsFactors = F
    )
    dir.tab(cbind(labs,pnames),paste0(filename,'Genes'),path)
  }
  if(plot.all) {
    xyarg <- append(xyarg,list(dat),0)
    col <- c('gray',col)
    pcex <- c(cex[1],pcex)
    ppch <- c(pch[1],ppch)
  }else legend <- legend[-1]
  xarg <- lapply(xyarg,'[',,'log2FoldChange.x')
  yarg <- lapply(xyarg,'[',,"log2FoldChange.y")
  args <- list(
    x=xarg,y=yarg,
    col = col,
    pch = ppch,cex=pcex
  )
  dir.eps(paste0(filename,'Gray'),path)
  do.call(plot,append(lapply(args,'[[',1),list(...)))
  dev.off()
  dir.eps(filename,path)
  mplot(
    args, ...
  )
  dev.off()
  dir.eps(paste0(filename,"Names"),path)
  mplot(args,legend = legend,...)
  do.call(text,append(pnames,list(offset=c(0.5,1),cex=.5)))
  dev.off()
  return(list(
    r=cor(xarg[[1]],yarg[[1]]),
    rho=cor(xarg[[1]],yarg[[1]],method = 'spearman')))
}

mplot <- function(
  args,fn=points,legend=NULL,legend.pos='topleft',
  xlim=range(args$x[[1]],na.rm = T),
  ylim=range(args$y[[1]],na.rm = T), 
  ...
){
  plot(NULL,xlim=xlim,ylim=ylim,...)
  if(!is.null(legend)) legend(legend.pos,legend=legend,col=args$col,pch=args$pch,cex=.75)
  # do.call(plot,append(lapply(args,'[[',1),list(...)))
  # args <- lapply(args,'[',-1)
  do.call(mapply,append(args,list(FUN=fn)))
  segments(xlim[1],0,xlim[2],0)
  segments(0,ylim[1],0,ylim[2])
}

