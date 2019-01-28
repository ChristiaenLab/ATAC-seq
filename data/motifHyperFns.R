combSetBinom <- function(motif,peakset,peaks,...){
  p.peakset <- length(intersect(peakset,peaks))/length(peaks)
  p.motif <- length(intersect(motif,peaks))/length(peaks)
  expected <- p.peakset*p.motif
  return(binom.test(
    length(Reduce(intersect,list(peakset,motif,peaks))),
    length(peaks),
    expected,
    alternative='greater'
  ))
  # p.sets <- mapply(
  #   function(x,y){
  #     length(intersect(sets[[x]],sets[[y]]))/length(sets[[y]])
  #   },
  #   1:(length(sets)-1),
  #   2:length(sets)
  # )
}

matchCombBinom <- function(matches,peaksets,...){
  peaks <- row.names(matches)
  motifs <- apply(matches,2,function(x) peaks[x])
  return(lapply(
    peaksets,
    function(peakset) lapply(
      motifs,combSetBinom,peakset=peakset,peaks=peaks,...
    )
  ))
}

# Accepts two logical matrices with experiments as columns and rows as observations.
# Columns in matches are vectors of which observations are successes.
# Columns in test are vectors used to select observations from matches to test enrichment
# runs hypergeometric test on each column of matches with parameters
# m=colSums(matches)
# n=nrow(matches)-m
# k=colSums(test)
# q=colSums(matches[test])
matchHyper <- function(matches,test,padj.method='fdr'){
  matches <- as.matrix(matches)
  peakct <- colSums(matches,na.rm = T)
  n.peakct <- nrow(matches)-peakct
  testct <- apply(test,2,function(x) colSums(matches[x,,drop=F],na.rm = T))
  k <- colSums(test,na.rm = T)
  # barplot(cbind(peakct/nrow(matches),testct/k))
  log2odds <- mapply(
    function(x,y) mapply(
      function(i,j) log2(
        (i/(y-i))/(j/(nrow(matches)-j))
      ),
      x,peakct
    ),
    as.data.frame(testct),k
  )
  testHyper <- mapply(function(q,k) mapply(
    phyper,q=q-1,k=k,m=peakct,n=n.peakct,lower.tail=F
  ),q=as.data.frame(testct),k=k)
  testFdr <- apply(testHyper,2,p.adjust,method=padj.method)
  colnames(log2odds) <- paste0('log2OddsRatio',colnames(testct))
  colnames(testHyper) <- paste0('p',colnames(testct))
  colnames(testFdr) <- paste0('fdr',colnames(testct))
  colnames(testct) <- paste0('counts',colnames(testct))
  res <- cbind(testct,log2odds,testHyper,testFdr)
  row.names(res) <- colnames(matches)
  return(cbind(counts=peakct,res))
}

matchFisher <- function(matches,test,padj.method='fdr'){
  matches <- as.matrix(matches)
  peakct <- colSums(matches,na.rm = T)
  n.peakct <- nrow(matches)-peakct
  testct <- apply(test,2,function(x) colSums(matches[x,],na.rm = T))
  k <- colSums(test,na.rm = T)
  # barplot(cbind(peakct/nrow(matches),testct/k))
  testHyper <- mapply(function(q,k) mapply(
    function(x,y,z,w) fisher.test(t(matrix(
      c(x,y,z-x,w-y),2,2
    )),alternative = "greater"),
    q,k-q,peakct,n.peakct,SIMPLIFY = F
  ), q=as.data.frame(testct),k=k,SIMPLIFY = F)
  p <- sapply(testHyper,sapply,'[[',"p.value")
  testFdr <- apply(p,2,p.adjust,method=padj.method)
  colnames(p) <- paste0('p',colnames(testct))
  odds <- log2(sapply(testHyper,sapply,'[[',"estimate"))
  
  res <- do.call(cbind,lapply(lapply(
    testHyper,sapply,"[",c("p.value","estimate")
  ),t))
  colnames(p) <- paste0('p',colnames(testct))
  colnames(testFdr) <- paste0('fdr',colnames(testct))
  colnames(testct) <- paste0('counts',colnames(testct))
  colnames(odds) <- paste0('log2OddsRatio',colnames(testct))
  res <- cbind(testct,odds,p,testFdr)
  row.names(res) <- colnames(matches)
  return(cbind(counts=peakct,res))
}

# runs matchHyper on the result of motifmatchr::matchMotifs
motifHyper <- function(matches,test,padj.method='fdr'){
  require(motifmatchr)
  tfs <- matches$name
  matches <- motifMatches(matches)
  colnames(matches) <- tfs
  # if(is.list(test)){
  #   test <- sapply(
  #     test, function(x) row.names(matches)%in%x
  #   )
  # }
  return(matchHyper(matches,test,padj.method))
}

# creates barplot of percent successes in each matches column recovered by each test column
# Bar groupings are matches columns. Colors are test columns.
# The first group of bars show percent of observations covered by each test set
# The first bar in each group shows percent of successes in matches.

# To compare treatments by feature, matches should be a matrix of significant peaks in each treatment.
# Test should be a matrix of which peaks are annotated to which features.
barplotHyper <- function(matches,test,p=T,...){
  peakct <- colSums(matches,na.rm = T)
  testct <- apply(test,2,function(x) colSums(matches[x,],na.rm = T))
  k <- colSums(test,na.rm = T)
  heights <- cbind(
    percent.peakome=c(peakome=NaN,k/nrow(matches)),
    rbind(peakome=peakct/nrow(matches),t(testct/peakct))
  )*100
  barplot(heights,beside = T,las=2,ylab = '% DA',legend.text = row.names(heights),ylim=c(0,max(heights[-1],na.rm=T)+4),...)
  if(p) {
    heights <- heights[-1,-1]
    fdr <- t(matchHyper(matches,test))
    fdr <- fdr[seq(2*nrow(heights)+2,3*nrow(heights)+1,1),]
    x <- sapply(
      seq(4.5+nrow(heights),by=nrow(heights)+2,length.out = ncol(heights)),
      function(x) seq(x,by=1,length.out = nrow(heights))
    )
    sig <- fdr<.05&fdr>=.01
    if(any(sig)) text(x[sig],heights[sig]+2,'*',srt=90)
    sig <- fdr<.01&fdr>=.001
    if(any(sig)) text(x[sig],heights[sig]+2,'**',srt=90)
    sig <- fdr<.001
    if(any(sig)) text(x[sig],heights[sig]+2,'***',srt=90)
  }
}

getGC <- function(peaks,genome=BSgenome.Cintestinalis.KH.KH2013,letters=c("C","G")){
  require(GenomicFeatures)
  require(BSgenome.Cintestinalis.KH.KH2013)
  if(class(peaks)=="GRanges") peaks <- split(peaks,seqnames(peaks))
  if(class(peaks)=="GRangesList"){
    strand(peaks)[strand(peaks)=="*"] <- "+"
    peaks <- extractTranscriptSeqs(genome,peaks)
  }
  peakct <- apply(letterFrequency(peaks,letters),2,sum)
  peakct <- peakct/sum(peakct)
  return(sum(peakct[c("C","G")]))
}

## A memory-efficient "letterFrequency" method for
## BSgenomeViews objects.
letterFrequency2 <- function(x, letters, OR="|",
                             as.prob=FALSE, ...)
{
  stopifnot(is(x, "BSgenomeViews"))
  chunksize <- 500000L
  chunks <- breakInChunks(length(x), chunksize)
  chunks <- as(chunks, "IRanges")
  ans_chunks <- lapply(seq_along(chunks),
                       function(i) {
                         x_chunk <- extractROWS(x, chunks[i])
                         letterFrequency(x_chunk, letters, OR=OR,
                                         as.prob=as.prob)
                       })
  do.call(rbind, ans_chunks)
}

getGC2 <- function(...){
  require(BSgenome.Cintestinalis.KH.KH2013)
  windowViews <- Views(BSgenome.Cintestinalis.KH.KH2013, windowRanges)
  gcFrequency <- letterFrequency2(windowViews, letters="GC", as.prob=TRUE)
}