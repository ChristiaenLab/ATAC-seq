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

reduceMotifs <- function(con,rmdup=T,khid.sub=T,selex.first=T){
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
  
  sel <- names(selex.pwm)%in%names(homer.pwm)
  names(selex.pwm)[sel] <- paste0(names(selex.pwm)[sel],"ANISEED")
  
  comb.pwm <- Reduce(append,list(selex.pwm,homer.pwm,cisbp.pwm))
  
  khToMotif <- rbind(
    cbind(ID=names(selex.pwm),GeneID=name(selex.pwm)),
    khToHomer[,c(1,3)],
    setNames(khToCisbp[,c(1,3)],c('ID',"GeneID")),
    stringsAsFactors=F
  )
  
  khToMotif <- khToMotif[khToMotif$ID%in%names(comb.pwm),]
  khToMotif <- khToMotif[!duplicated(khToMotif),]
  
  if(selex.first){
    khToMotif <- khToMotif[!khToMotif$GeneID%in%name(selex.pwm),]
    motifs <- selex.pwm
    names(motifs) <- name(selex.pwm)
  } else motifs <- NULL
  # khToHomer <- split(khToHomer$GeneID,khToHomer$ID)
  # khToHomer <- khToHomer[sapply(khToHomer,function(x) !any(x%in%name(selex.pwm)))]
  # 
  # homer.pwm <- homer.pwm[names(khToHomer)]
  # names(homer.pwm) <- sapply(khToHomer,paste,collapse=';')
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
  
  
  motifs <- Reduce(append,list(motifs,tmp2))
  if(rmdup){
    motifs <- motifs[!duplicated(names(motifs))]
  }
  ann <- getAnnotation(con)
  motifs <- nameMotifs(motifs,ann$gene.names,khid.sub=khid.sub)
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

lowerGeneName <- function(x){
  if(
    grepl("[A-Z]{2}",x)&!(grepl("^KH\\.",x)|grepl("^SI:",x))
  ){
    tmp <- sub(".*?([A-Z])([A-Z]+)",'\\2',x)
    return(sub(tmp,tolower(tmp),x))
  } else {x}
}

khToName <- function(x,gene.names){
  x <- gene.names[x,"UniqueNAME"]
  x <- sub(
    '^KH\\.[A-Z][0-9]+\\.[0-9]+_','',sub("KH2013:",'',x)
  )
  x <- sub(';$','',x)
  x <- sapply(x,lowerGeneName)
  return(x)
}

nameMotifs <- function(motifs,gene.names,khid.sub=T){
  require(TFBSTools)
  tf.family <- sapply(tags(motifs),'[[',"Family_Name")
  tf.family[grep("(Gata|Zn?F)",tf.family,T)] <- "Zinc finger"
  tf.family[tf.family=="NR"] <- "Nuclear receptor"
  tf.family[tf.family=="ETS"] <- "Ets"
  tf.family[tf.family=="MAD"] <- "SMAD"
  tf.family[tf.family=="MADS"] <- "MADS-box"
  tf.family[tf.family=="Paired"] <- "Homeodomain,POU"
  tf.family[tf.family=="Homeodomain,Paired box"] <- "Paired box"
  tf.family[tf.family=="Paired,Homeobox"] <- "Paired box"
  tf.family[tf.family=="EBF"] <- "HLH"
  tf.family[tf.family=="POU,Homeobox,HMG"] <- "HMG"
  tf.family[tf.family=="CTF,Forkhead"] <- "CTF"
  tf.family[tf.family=="ETS:IRF"] <- "IRF"
  tf.family[tf.family=="promoter"] <- "TBP"
  tf.family[tf.family=="Homeobox,bHLH"] <- "Homeodomain"
  tf.family[tf.family=="AP2"] <- "AP-2"
  tf.family[tf.family=="E2F/TDP"] <- "E2F"
  
  tf.tags <- mapply(function(x,y){
    x$Family_Name <- y
    return(x)
  },tags(motifs),tf.family)
  
  tf.name <- sapply(tags(motifs),'[[',"DBID.1")
  tf.kh <- names(motifs)
  
  tf.kh.gene <- strsplit(tf.kh,';')
  tf.kh.gene <- lapply(tf.kh.gene,khToName,gene.names)
  tf.kh.gene <- lapply(tf.kh.gene,function(x) x[!duplicated(x)])
  
  if(khid.sub){
    sel <- sapply(tf.kh.gene,function(x) all(grepl("^KH\\.[A-Z][0-9]+",x)))
    tf.kh.gene[sel] <- sapply(ID(motifs)[sel],list)
  }
  
  tf.kh.gene <- lapply(tf.kh.gene,sub,pattern='V\\$',replacement='')
  
  # tf.kh.gene <- lapply(tf.kh.gene,function(y) sapply(y,lowerGeneName))#function(x) if(
  #   grepl("^[A-Z]{2}",x)&!(grepl("^KH\\.",x)|grepl("^SI:",x))
  # ){
  #   tmp <- sub("(^.)([A-Z]+)",'\\2',x)
  #   return(sub(tmp,tolower(tmp),x))
  # } else {x}))
  
  tf.kh.gene <- lapply(tf.kh.gene,Reduce,f=mergeGeneName)
  
  tf.kh.gene <- sapply(tf.kh.gene,paste,collapse=';')
  
  tf.family[tf.kh.gene=="Ctcf"] <- "Zinc finger"
  tf.family[tf.kh.gene=="Rbpj"] <- "CSL"
  
  motifs <- mapply(
    PWMatrix,
    names(motifs),
    tf.kh.gene,
    strand="*",
    tags=tf.tags,
    profileMatrix=Matrix(motifs)
  )
  names(motifs) <- make.unique(tf.kh.gene)
  
  motifs <- do.call(PWMatrixList,motifs)
  return(motifs)
}