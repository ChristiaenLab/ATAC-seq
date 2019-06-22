getMotifs <- function(con){
  require(TFBSTools)
  require(motifmatchr)
  require(DBI)
  require(BSgenome.Cintestinalis.KH.JoinedScaffold)
  
  source('data/chromVarFns.R')
  source('data/sqlfns.R')
  source('data/getSelex.R')
  
  khToHomer <- dbReadTable(con,'homer_orthologs')
  khToCisbp <- dbReadTable(con,'cisbp_orthologs')
  
  selex.pwm <- getSelex()
  cisbp.pwm <- getCisbpMotifs()
  homer.pwm <- getHomerMotifs('known.motifs')
  
  khToHomer <- split(khToHomer$GeneID,khToHomer$ID)
  khToHomer <- khToHomer[sapply(khToHomer,function(x) !any(x%in%name(selex.pwm)))]
  
  homer.pwm <- homer.pwm[names(khToHomer)]
  names(homer.pwm) <- sapply(khToHomer,paste,collapse=';')
  names(selex.pwm) <- name(selex.pwm)
  
  khToCisbp.selex <- khToCisbp[khToCisbp$Motif_Type=="SELEX",]
  khToCisbp.selex <- split(khToCisbp.selex$GeneID,khToCisbp.selex$Motif_ID)
  khToCisbp.selex <- khToCisbp.selex[sapply(
    khToCisbp.selex,
    function(x) !any(x%in%c(names(selex.pwm),unlist(khToHomer)))
  )]
  cisbp.selex.pwm <- cisbp.pwm[names(khToCisbp.selex)]
  names(cisbp.selex.pwm) <- sapply(khToCisbp.selex,paste,collapse=';')
  
  khToCisbp.other <- split(khToCisbp$GeneID,khToCisbp$Motif_ID)
  khToCisbp.other <- khToCisbp.other[names(khToCisbp.other)%in%names(cisbp.pwm)]
  khToCisbp.other <- khToCisbp.other[sapply(
    khToCisbp.other,
    function(x) !any(x%in%c(names(selex.pwm),unlist(khToHomer),unlist(khToCisbp.selex)))
  )]
  cisbp.pwm <- cisbp.pwm[names(khToCisbp.other)]
  names(cisbp.pwm) <- sapply(khToCisbp.other,paste,collapse=';')
  
  motifs <- Reduce(append,list(selex.pwm,homer.pwm,cisbp.selex.pwm,cisbp.pwm))
  return(motifs)
}

reduceMotifs <- function(con){
  require(TFBSTools)
  require(motifmatchr)
  require(DBI)
  require(BSgenome.Cintestinalis.KH.JoinedScaffold)
  
  source('data/chromVarFns.R')
  source('data/sqlfns.R')
  
  khToHomer <- dbReadTable(con,'homer_orthologs')
  khToCisbp <- dbReadTable(con,'cisbp_orthologs')
  khToMotif <- rbind(khToHomer[,c(1,3)],setNames(khToCisbp[,c(1,3)],c('ID',"GeneID")))
  
  selex.pwm <- getSelex()
  cisbp.pwm <- getCisbpMotifs()
  homer.pwm <- getHomerMotifs('known.motifs')
  comb.pwm <- append(homer.pwm,cisbp.pwm)
  khToMotif <- khToMotif[khToMotif$ID%in%names(comb.pwm),]
  khToMotif <- khToMotif[!khToMotif$GeneID%in%name(selex.pwm),]
  khToMotif <- khToMotif[!duplicated(khToMotif),]
  
  # khToHomer <- split(khToHomer$GeneID,khToHomer$ID)
  # khToHomer <- khToHomer[sapply(khToHomer,function(x) !any(x%in%name(selex.pwm)))]
  # 
  # homer.pwm <- homer.pwm[names(khToHomer)]
  # names(homer.pwm) <- sapply(khToHomer,paste,collapse=';')
  names(selex.pwm) <- name(selex.pwm)
  # 
  # khToCisbp <- khToCisbp[!khToCisbp$GeneID%in%c(
  #   names(selex.pwm),unlist(khToHomer)
  # )&khToCisbp$Motif_ID%in%names(cisbp.pwm),]
  fn <- function(x,y,res=NULL){
    motifs <- split(x$khToMotif$GeneID,x$khToMotif$ID)
    tmp <- motifs[sapply(motifs,length)<=y]
    if(length(tmp)>0){
      res <- append(res,tmp)
      x$khToMotif <- x$khToMotif[!x$khToMotif$GeneID%in%unlist(res),]
      fn(x,y,res)
    }else{
      x$motifs <- append(x$motifs,res)
      return(x)
    }
  }
  tmp <- Reduce(
    fn,
    1:max(sapply(split(khToMotif$GeneID,khToMotif$ID),length)),
    list(motifs=list(),khToMotif=khToMotif)
  )
  tmp2 <- comb.pwm[names(tmp$motifs)]
  names(tmp2) <- sapply(tmp$motifs,paste,collapse=';')
  # tmp2 <- Reduce(function(x,y){
  #   motifs <- comb.pwm[names(y)]
  #   names(motifs) <- sapply(y,paste,collapse=';')
  #   return(append(x,motifs))
  # },tmp$motifs,NULL)
  
  
  motifs <- Reduce(append,list(selex.pwm,tmp2))
  return(motifs)
}

mergeGeneName <- function(x,y){
  expr <- sub('^([A-Za-z][A-Za-z0-9]*[A-Za-z])[0-9/]+$','\\1',x)
  sel <- sapply(paste0("^",expr,"[0-9/]+$"),grepl,y)
  if(any(sel)){
    nums <- as.numeric(unlist(strsplit(sub(expr,'',x),'/')))
    numsToAdd <- as.numeric(unlist(strsplit(sub(expr,'',y),'/')))
    nums <- unique(c(nums,numsToAdd))
    nums <- nums[order(nums)]
    y <- paste(as.character(nums),collapse='/')
    x[sel] <- paste0(expr,y)
  }else{
    x <- append(x,y)
  }
  return(x)
}

nameMotifs <- function(motifs,gene.names){
  require(TFBSTools)
  tf.name <- sapply(tags(motifs),'[[',"DBID.1")
  tf.kh <- names(motifs)
  
  tf.kh.gene <- strsplit(tf.kh,';')
  tf.kh.gene <- lapply(tf.kh.gene,function(x) sub(
    '^KH\\.[A-Z][0-9]+\\.[0-9]+_','',sub("KH2013:",'',gene.names[x,'UniqueNAME'])
  ))
  tf.kh.gene <- sapply(tf.kh.gene,sub,pattern=';$',replace='')
  tf.kh.gene <- lapply(tf.kh.gene,function(x) x[!duplicated(x)])
  sel <- sapply(tf.kh.gene,function(x) all(grepl("^KH\\.[A-Z][0-9]+",x)))
  tf.kh.gene[sel] <- sapply(ID(motifs)[sel],list)
  tf.kh.gene <- lapply(tf.kh.gene,sub,pattern='V\\$',replacement='')
  tf.kh.gene <- lapply(tf.kh.gene,function(y) sapply(y,function(x) if(
    grepl("^[A-Z]{2}",x)&!(grepl("^KH\\.",x)|grepl("^SI:",x))
  ){
    tmp <- sub("(^.)([A-Z]+)",'\\2',x)
    return(sub(tmp,tolower(tmp),x))
  } else {x}))
  
  tf.kh.gene <- lapply(tf.kh.gene,Reduce,f=mergeGeneName)
  
  tf.kh.gene <- sapply(tf.kh.gene,paste,collapse=';')
  names(motifs) <- make.unique(tf.kh.gene)
  sel <- grep("^KH\\.[A-Z][0-9]+",names(motifs))
  names(motifs)[sel] <- tf.name[sel]
  return(motifs)
}