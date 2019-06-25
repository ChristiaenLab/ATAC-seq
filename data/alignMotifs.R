alignMotifs <- function(fasta,motifs,species=c("Crobusta","Csavignyi"),suffix='',family="Family_Name",...){
  require(rtracklayer)
  require(motifmatchr)
  require(TFBSTools)
  dat <- import(fasta,'fasta')
  names(dat) <- species
  family <- unlist(lapply(tags(motifs),'[[',family))
  matches <- matchMotifs(
    motifs,
    dat,
    out='positions',
    ...
  )
  matchpos <- sapply(1:length(species),function(x){
    cr <- lapply(matches,'[[',x)
    sel <- sapply(cr,length)>0
    family <- family[sel]
    cr <- cr[sel]
    cr <- mapply(function(x,y){
      mcols(x)$family <- y
      return(x)
    },cr,family)
    cr <- do.call(IRangesList,cr)
    cr <- unlist(cr)
    cr <- GRanges(
      species[x],
      do.call(IRanges,as.data.frame(cr)),
      mcols(cr)$strand,
      score=mcols(cr)$score,
      family=mcols(cr)$family
    )
    cr$seq <- as.character(Views(dat[[x]],ranges(cr)))
    return(cr)
  })
  names(matchpos) <- species
  matches <- Reduce(c,matchpos)
  # names(matches) <- name(motifs[names(matches)])
  file <- paste0(sub('\\..*$','',fasta),suffix)
  dir.export(matches,paste0(file,"Matches"),format = "gff3")
  
  x <- sapply(matchpos,function(x) sapply(
    split(x,names(x)),
    function(y) mcols(y)[which.max(y$score),c('score','family')]
  ))
  x <- lapply(x,do.call,what=rbind)
  
  x <- Reduce(function(i,j) merge(i,j,c("row.names",'family'),all=T),x)
  y <- setNames(x[,-1:-2],species)
  y <- as.matrix(y)
  row.names(y) <- x[,1]#sub('(.{,20}).*','\\1',x[,1])
  # y[is.na(y)] <- 0
  sel <- !apply(y,1,anyNA)
  plotPeakMatches(
    t(y[sel,]),
    paste0(file,'Conserved'),
    x$family[sel]
  )
}
