GRexpand <- function(x,start,end=start){
  x@ranges <- IRanges(
    x@ranges@start-start,
    width=x@ranges@width+start+end)
  return(x)
}

bedGR <- function(x,...) {
  x <- as.list(read.table(x,...))
  x$seqnames <- Rle(x[[1]])
  x$ranges <- IRanges(as.numeric(x[[2]]),as.numeric(x[[3]]))
  x$strand <- Rle(rep('*',length(x[[3]])))
  return(do.call(GRanges,x))
}

GRbed <- function(x,filename,path){
  x <- cbind(
    as.vector(seqnames(x)),
    ranges(x)@start,ranges(x)@start+ranges(x)@width,
    as.vector(strand(x)),as.data.frame(mcols(x)))
  sapply(2:3,function(y) x[x[,y]<0,y] <<- 0)
  x <- x[!apply(
    x[,2:3]==0,1,
    function(y) do.call('&',as.list(y))
    ),] 
  path <- paste(path, Sys.Date(), sep = '/')
  if(!dir.exists(path)) dir.create(path,recursive = T)
  filename <- paste(path, filename, sep = '/')
  write.table(x,paste(filename,'bed',sep='.'),sep='\t', quote = F,col.names = F,row.names = F)
}

readPeakGenes <- function(file){
  peaks <- read.table(file,colClasses = 'character')
  peaks$peak <- apply(peaks[,1:3],1,paste,collapse=':')
  peaks <- sapply(
    unique(peaks$peak),function(x) paste(unique(peaks[peaks$peak==x,9]),collapse = ';'))
  peaks <- cbind(peaks,sapply(
    peaks,
    function(x) paste0(
      gene.names[unlist(strsplit(as.character(x),';')),],collapse = ';')))
  return(peaks)
}

readGTF <- function(gtf){
  gtf <- read.delim(gtf, header=FALSE)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",     
                     "attributes")
  gtf$GeneID <- as.factor(sub(".*gene_id\\s*(.*);", "\\1", gtf$attributes))
  # gtf$featureID <- paste0(gtf$GeneID,';',gtf$feature)
  # gtf <- gtf[order(gtf$seqname, gtf$start), ]
  # gtf1 <- gtf1[!duplicated(gtf1$GeneID), ] # Returns only rows with most left gene positions.
  # gtf2 <- gtf[order(gtf$seqname, -gtf$start), ]
  # gtf2 <- gtf2[!duplicated(gtf2$GeneID), ] # Returns only rows with most right gene positions.
  # gtf1 <- gtf1[order(gtf1$GeneID), ]
  # gtf2 <- gtf2[order(gtf2$GeneID), ]
  # gtf1[,5] <- gtf2[,5] # Assigns proper end postions
  # gtfgenepos <- gtf1 # Returns gtf containing only start/end positions for all genes
  # gtfgenepos <- gtfgenepos[order(gtfgenepos$seqname, gtfgenegtfpos$start), ]
  gtf <- makeGRangesListFromDataFrame(gtf,'GeneID',keep.extra.columns=T)
  
  # cds <- gtf[]
  # cds <- reduce(gtf[grep('CDS',names(gtf))])
  # names(cds) <- sub(';CDS','',names(cds))
  # utr <- range(gtf[grep('exon',names(gtf))])
  # names(utr) <- sub(";exon",'',names(utr))
  # introns <- disjoin(cds)
  # utr <- 
  return(gtf)
}

subsetGRList <- function(x,expr){
  sel <- eval(expr,elementMetadata(unlist(x)))
  sel <- relist(sel,x)
  res <- lapply(1:length(x), function(i) x[[i]][sel[[i]]])
  res <- do.call(GRangesList,res)
  names(res) <- names(x)
  return(res)
}

parseGTF <- function(gtf){
  require(GenomicRanges)
  gtf <- readGTF(gtf)
  
  gtf <- reduce(gtf)
  # cds <- subsetGRList(gtf,expression(feature=='CDS'))
  # utr <- subsetGRList(gtf,expression(feature=='exon'))
  introns <- disjoin(gtf)
  genes <- unlist(range(gtf))
  intergenic <- gaps(genes)
  
  # gtfgaps <- gaps(gtf)
  # genes <- sapply(unique(gtf@elementMetadata$GeneID),function(x) range(gtf[gtf@elementMetadata$GeneID==x]))
}

getFeatures <- function(
  gff,
  tssflank=c(107,107),
  ttsflank=c(200,200),
  promoterflank=c(promoter500=500,promoter1k=500),
  tsswindow=10000,ttswindow=10000
){
  require(GenomicRanges)
  require(GenomicFeatures)
  gff <- gff[sapply(elementMetadata(gff)$Parent,function(x) length(x)>0)]
  elementMetadata(gff)$GeneID <- unlist(sub('\\.v.*','',elementMetadata(gff)$Parent))
  elementMetadata(gff)$GeneName <- gene.names[elementMetadata(gff)$GeneID,]
  txdb <- makeTxDbFromGRanges(gff)
  
  # tpt <- transcriptsBy(txdb)
  # fputr <- fiveUTRsByTranscript(txdb)
  # tputr <- threeUTRsByTranscript(txdb)
  # cds <- cdsBy(txdb,'gene')
  # exons <- exonsBy(txdb,'gene')
  introns <- reduce(unlist(intronsByTranscript(txdb)))
  
  fputr <- gff[elementMetadata(gff)$type=='five_prime_UTR']
  fputr <- split(fputr,elementMetadata(fputr)$GeneID)
  tputr <- gff[elementMetadata(gff)$type=='three_prime_UTR']
  tputr <- split(tputr,elementMetadata(tputr)$GeneID)
  cds <- gff[elementMetadata(gff)$type=='CDS']
  cds <- split(cds,elementMetadata(cds)$GeneID)
  gene <- gff[elementMetadata(gff)$type%in%c("CDS",'five_prime_UTR','three_prime_UTR')]
  gene <- split(gene,elementMetadata(gene)$GeneID)
  genebody <- unlist(range(gene,ignore.strand=T))
  # geneGaps <- gaps(reduce(unlist(gene),ignore.strand=T))
  # introns <- first(findOverlapPairs(geneGaps,genebody,type='within'))
  # introns <- mapply(gaps,reduce(gene,ignore.strand=T),start=start(genebody),end=end(genebody))
  # introns <- do.call(GRangesList,introns)
  
  tpt <- unlist(range(split(gff,elementMetadata(gff)$ID)))
  elementMetadata(tpt)$GeneID <- sub('\\.v.*','',names(tpt))
  tpt <- split(tpt,elementMetadata(tpt)$GeneID)
  promoterGene <- reduce(promoters(
    tpt,
    tssflank[1]+sum(promoterflank),
    tssflank[2]
  ))
  
  # cds <- reduce(cdsBy(txdb,'gene'))
  # exon <- exonsBy(gff,'gene')
  # tpt <- transcriptsBy(txdb,'gene')
  # intergenic <- gaps(tpt)
  tss <- promoters(txdb,tssflank[1],tssflank[2])
  dir.export(promoters(
    txdb,
    tssflank[1]+sum(promoterflank),
    tssflank[2]
  ),'promoters',format = 'gtf')
  elementMetadata(tss)$GeneID <- sub('\\.v.*','',elementMetadata(tss)$tx_name)
  dir.export(tss,'tssCenter',path)
  tss <- split(tss,elementMetadata(tss)$GeneID)
  tss <- reduce(tss)
  dir.export(tss,'TSS',path)
  
  tts <- flank(tpt,ttsflank[1],start = F)
  tts <- resize(tts,sum(ttsflank),fix = 'end')
  
  tsswindow <- flank(unlist(range(tpt)),tsswindow)
  ttswindow <- flank(unlist(range(tpt)),ttswindow,start = F)
  promoters <- setNames(
    Reduce(flank,promoterflank,accumulate = T,init = tpt)[-1],
    names(promoterflank)
  )
  promoters <- promoters[length(promoters):1]
  # promoter500 <- flank(tpt,500)
  # promoter1k <- flank(promoter500,500)
  intergenic <- gaps(reduce(union(
    unlist(Reduce(
      union,
      append(list(tts,tss),promoters)
    )),
    genebody),ignore.strand=T))
  
  return(list(
    genebody=genebody,
    promoterGene=promoterGene,
    features=append(
      promoters,
      list(
        TSS=tss,five_prime_utr=fputr,
        CDS=cds,intron=introns,TTS=tts,
        three_prime_utr=tputr,intergenic=intergenic
      )
    )
  ))
}

getAnn <- function(query,subject,colname){
  ann <- findOverlaps(query,subject)
  ann <- split(ann@to,ann@from)
  ann <- sapply(ann,function(x) unique(names(subject)[x]))
  ann <- sapply(ann,paste,collapse = ';')
  elementMetadata(query)[,colname] <- NA
  elementMetadata(query)[,colname][as.numeric(names(ann))] <- unlist(ann)
  return(query)
}

annotatePeaks <- function(peaks,genes,window=10000,features=NULL){
  require(GenomicRanges)
  genes <- resize(genes,width(genes)+2*window,'center')
  annotation <- findOverlaps(peaks,genes)
  # names(peaks) <- paste(seqnames(peaks),start(peaks)-1,end(peaks),sep = ':')
  # names(peaks) <- Reduce(
  #   function(x,y){
  #     x <- unlist(strsplit(x,'\\.'))
  #     if(y==x[1]) {
  #       return(paste0(y,'.',as.character(as.numeric(x[2])+1)))
  #     }else return(paste0(y,'.',as.character(1)))
  #   },
  #   seqnames(peaks),'KhC1.0',
  #   accumulate=T
  # )[-1]
  # seqnames(peaks) <- droplevels(seqnames(peaks))
  names(peaks) <- unlist(mapply(
    paste0,
    as.character(levels(droplevels(seqnames(peaks)))),'.',
    mapply(seq,1,table(droplevels(seqnames(peaks))))
  ))
  res <- list(
    annotation=annotation,
    peaks=peaks,
    genes=genes,
    window=window
  )
  res$peakToGene <- GRpeakToGene(res)
  res$geneToPeak <- GRgeneToPeak(res)
  if(!is.null(features)){
    res$features <- sapply(
      features,
      function(x) overlapsAny(peaks,unlist(x))
    )
    res$features <- cbind(
      data.frame(
        chr=seqnames(peaks),
        start=start(peaks),
        end=end(peaks),
        row.names = names(peaks)
      ),res$features
    )
    # row.names(res$features) <- names(peaks)
  }
  return(res)
}

splitBy <- function(from,to,names,elements){
  res <- split(elements[to],names[from])
  res <- sapply(res,unique)
  # res <- sapply(elements,function(x) elements[x])
  # names(res) <- names[as.numeric(names(res))]
  return(res)
}

GRpeakToGene <- function(ann,peaks=names(ann$peaks)){
  genes <- split(ann$annotation@to,ann$annotation@from)
  genes <- sapply(genes,unique)
  genes <- sapply(genes,function(x) names(ann$genes)[x])
  names(genes) <- names(ann$peaks)[as.numeric(names(genes))]
  return(genes)
}

GRgeneToPeak <- function(ann,peaks=names(ann$peaks)){
  peaks <- split(ann$annotation@from,ann$annotation@to)
  peaks <- sapply(peaks,unique)
  peaks <- sapply(peaks,function(x) names(ann$peaks)[x])
  names(peaks) <- names(ann$genes)[as.numeric(names(peaks))]
  return(peaks)
}

# runs a binomial test for enrichment of peaks in annotation associated to genes
# compared to bg
# the expected probability is the fraction of bg covered by genes
# genes should be a vector of GeneIDs
# annotation should be a list of PeakID vectors with the GeneID as names(annotation)
# bg should be a GRanges of genes
geneSizeBinom <- function(
  genes,annotation,bg,
  alternative='two.sided',conf.level=.99,...
){
  windowSize <- sum(width(reduce(bg)))
  query <- sapply(genes,intersect,names(bg),simplify = F)
  querySize <- sapply(query,function(i) sum(width(reduce(bg[i]))))
  x <- sapply(query,function(i) length(unique(unlist(peakGeneAnnotation$geneToPeak[i]))))
  n <- length(unique(unlist(peakGeneAnnotation$geneToPeak[names(bg)])))
  p <- querySize/windowSize
  res <- mapply(
    binom.test,x,n,p,
    alternative=alternative,conf.level=conf.level,...,
    SIMPLIFY = F
  )
  estimate <- sapply(res,'[[','estimate')
  lowerConf <- sapply(res,'[[','conf.int')[1,]
  upperConf <- sapply(res,'[[','conf.int')[2,]
  dat <- cbind(
    n=x,
    p=sapply(res,'[[','p.value'),
    padj=p.adjust(sapply(res,'[[','p.value')),
    null.value=p,
    estimate=estimate,
    lowerConf=lowerConf,
    upperConf=upperConf,
    log2oddsRatio=log2(estimate/p),
    log2oddsLower=log2(lowerConf/p),
    log2oddsUpper=log2(upperConf/p)
  )
  row.names(dat) <- names(res)
  return(dat)
}

# getIntrons <- function(genes) sapply(genes,function(x) gaps(x,start = start(x)[1]))