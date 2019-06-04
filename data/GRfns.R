getFeatures <- function(
# accepts a GFF file name
# returns a list of GRangesLists corresponding to genomic features
  gff,
  gene.names,
  tssflank=c(107,107),
  ttsflank=c(200,200),
  promoterflank=c(promoter500=500,promoter1k=500),
  tsswindow=10000,ttswindow=10000
){
  require(GenomicRanges)
  require(GenomicFeatures)
  gff <- import(gff)
  gff <- gff[sapply(elementMetadata(gff)$Parent,function(x) length(x)>0)]
  elementMetadata(gff)$GeneID <- unlist(sub('\\.v.*','',elementMetadata(gff)$Parent))
  elementMetadata(gff)$GeneName <- gene.names[elementMetadata(gff)$GeneID,]
  txdb <- makeTxDbFromGRanges(gff)
  
  introns <- intronsByTranscript(txdb)
  
  fputr <- gff[elementMetadata(gff)$type=='five_prime_UTR']
  fputr <- split(fputr,elementMetadata(fputr)$GeneID)
  tputr <- gff[elementMetadata(gff)$type=='three_prime_UTR']
  tputr <- split(tputr,elementMetadata(tputr)$GeneID)
  cds <- gff[elementMetadata(gff)$type=='CDS']
  cds <- split(cds,elementMetadata(cds)$GeneID)
  gene <- gff[elementMetadata(gff)$type%in%c("CDS",'five_prime_UTR','three_prime_UTR')]
  gene <- split(gene,elementMetadata(gene)$GeneID)
  genebody <- unlist(range(gene,ignore.strand=T))
  
  tpt <- unlist(range(split(gff,elementMetadata(gff)$ID)))
  elementMetadata(tpt)$GeneID <- sub('\\.v.*','',names(tpt))
  tpt <- split(tpt,elementMetadata(tpt)$GeneID)
  promoterGene <- reduce(promoters(
    tpt,
    tssflank[1]+sum(promoterflank),
    tssflank[2]
  ))
  
  tss <- promoters(txdb,tssflank[1],tssflank[2])
  dir.export(promoters(
    txdb,
    tssflank[1]+sum(promoterflank),
    tssflank[2]
  ),'promoters',format = 'gtf')
  elementMetadata(tss)$GeneID <- sub('\\.v.*','',elementMetadata(tss)$tx_name)
  dir.export(tss,'tssCenter',format = 'gtf')
  tss <- split(tss,elementMetadata(tss)$GeneID)
  tss <- reduce(tss)
  dir.export(tss,'TSS')
  
  tts <- flank(tpt,ttsflank[1],start = F)
  tts <- resize(tts,sum(ttsflank),fix = 'end')
  
  tsswindow <- flank(unlist(range(tpt)),tsswindow)
  ttswindow <- flank(unlist(range(tpt)),ttswindow,start = F)
  promoters <- setNames(
    Reduce(flank,promoterflank,accumulate = T,init = tpt)[-1],
    names(promoterflank)
  )
  promoters <- promoters[length(promoters):1]
  intergenic <- GRangesList(gaps(reduce(union(
    unlist(Reduce(
      union,
      append(list(tts,tss),promoters)
    )),
    genebody),ignore.strand=T)))
  
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

annotatePeaks <- function(peaks,genes,window=10000,features=NULL){
  # peaks is a GRanges of peaks 
  # genes is a GRanges of gene loci
  # window is the number of bps genes will be expanded in either direction
  # features is a named list of GRangesLists
  
  # returns a list with attributes:
    # annotation  Hits object for overlaps between peaks and genes
    # peaks   the input peaks
    # genes   the input genes
    # window  the input genes expanded by window
    # peakToGene  a list of GeneID vectors split by PeakID
    # geneToPeak  a list of PeakID vectors split by GeneID
    # features  a data. frame with columns indicating whether each peak overlaps with each feature in features
  require(GenomicRanges)
  genes <- resize(genes,width(genes)+2*window,'center')
  annotation <- findOverlaps(peaks,genes)
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
      overlapsAny,query=peaks
    )
    res$features <- cbind(
      data.frame(
        chr=seqnames(peaks),
        start=start(peaks),
        end=end(peaks),
        row.names = names(peaks)
      ),res$features
    )
  }
  return(res)
}

splitBy <- function(from,to,names,elements){
  # wrapper function for split ensuring each vector contains only unique values
  res <- split(elements[to],names[from])
  res <- sapply(res,unique)
  return(res)
}

GRpeakToGene <- function(ann,peaks=names(ann$peaks)){
  # splits a findOverlaps result with peaks as the query and genes as the subject by peaks 
  genes <- split(ann$annotation@to,ann$annotation@from)
  genes <- sapply(genes,unique)
  genes <- sapply(genes,function(x) names(ann$genes)[x])
  names(genes) <- names(ann$peaks)[as.numeric(names(genes))]
  return(genes)
}

GRgeneToPeak <- function(ann,peaks=names(ann$peaks)){
  # splits a findOverlaps result with peaks as the query and genes as the subject by genes
  peaks <- split(ann$annotation@from,ann$annotation@to)
  peaks <- sapply(peaks,unique)
  peaks <- sapply(peaks,function(x) names(ann$peaks)[x])
  names(peaks) <- names(ann$genes)[as.numeric(names(peaks))]
  return(peaks)
}

geneSizeBinom <- function(
  genes,annotation,bg,
  alternative='two.sided',conf.level=.99,...
){
# runs a binomial test for enrichment of peaks in annotation associated to genes
# compared to bg
# the expected probability is the fraction of bg covered by genes
# genes should be a vector of GeneIDs
# annotation should be a list of PeakID vectors with the GeneID as names(annotation)
# bg should be a GRanges of genes
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