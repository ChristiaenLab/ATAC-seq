addSubdir <- function(dir,subdir){
  dir <- paste0(dir,'/',subdir,'/')
  if(!dir.exists(dir))dir.create(dir,recursive = T)
  return(dir)
}

pathDate <- function(
  date=Sys.Date(),out='CISBP/Ciona_intestinalis_2017_11_01_3_41_pm/chromVAR/'
) sapply(
  list.dirs(out,recursive=F),
  paste0,'/',date
)

getHomerMotifs <- function(homerdat){
  motifs <- readLines(homerdat)
  motifs <- split(
    motifs,
    Reduce(function(x,y) if(y) c(x,x[length(x)]+1) else c(x,x[length(x)]),grepl('^>',motifs))
  )
  require(TFBSTools)
  motifs <- lapply(motifs,strsplit,'\t')
  homer.motifs <- do.call(PWMatrixList,lapply(
    motifs,
    function(x){
      profileMatrix <- sapply(x[-1],as.numeric)
      profileMatrix <- t(t(profileMatrix)/apply(profileMatrix,2,sum))
      row.names(profileMatrix) <- c("A","C","G","T")
      # ID <- sub('^>','',x[[1]][1])
      tags <- strsplit(x[[1]][2],'/')[[1]]
      tags <- c(strsplit(tags[1],'[()]')[[1]],tags[-1])
      DBID.1 <- tags[1]
      ID <- make.names(DBID.1,T)
      Family_Name <- if(
        is.na(tags[2])|grepl("\\?",tags[2])|grepl("\\.",tags[2])
      ){"NA"} else if(
        grepl("Bias",tags[2],ignore.case = T)|grepl("repeat",tags[2],ignore.case = T)
      ){"SeqBias"} else if(
        grepl("promoter",tags[2],ignore.case = T)
      ){"promoter"} else tags[2]
      return(PWMatrix(
        ID=ID,
        name=ID,
        profileMatrix = profileMatrix,
        tags = list(Family_Name=Family_Name,DBID.1=DBID.1)
      ))
    }
  ))
  names(homer.motifs) <- make.names(ID(homer.motifs),T)
  return(homer.motifs)
}

getCisbpMotifs <- function(cisbp,promoters=c(
  "BREd","BREu","DCEI-DCEIII","DPE","DRE","E-box","Inr_fly","Inr_human",
  "ohler","Pause_button",'TATA-box',"TCT","XCPE","MTE"
)){
  require(TFBSTools)
  motifs <- lapply(
    list.files(
      paste0('CISBP/',cisbp,'/pwms_all_motifs/'),
      full.names = T),
    function(x){
      mat <- t(as.matrix(read.delim(x)[,c("A","C","G","T")]))
      if(dim(mat)[2]>0){
        sapply(1:ncol(mat),function(y){
          mat[which.max(mat[,y]),y] <<- mat[which.max(mat[,y]),y]-(sum(mat[,y])-1)
        })
        # mat <- mat/apply(mat,2,sum)
        return(list(ID=sub('.txt','',sub('.*\\/','',x)),profileMatrix = mat,strand='*'))
      }
    })
  motifs <- motifs[sapply(motifs,class)=="list"]
  
  motifid <- unlist(sapply(motifs,'[','ID'))
  
  motifdat <- read.delim(
    paste0('CISBP/',cisbp,'/TF_Information_all_motifs_plus.txt'),
    stringsAsFactors = F)
  addmotif <- motifid[!motifid%in%motifdat$Motif_ID]
  sapply(addmotif, function(x) {
    motifdat[x,c('Motif_ID',"DBID.1",'Family_Name')] <<- x
    motifdat[x,'Family_Name'] <<- sub('[0-9]+$','',motifdat[x,'Family_Name'])
    if(
      motifdat[x,'Family_Name']%in%promoters
    ) motifdat[x,'Family_Name'] <<- "promoter"
  })
  motifdat <- motifdat[!duplicated(motifdat$Motif_ID),]
  row.names(motifdat) <- motifdat$Motif_ID
  motifdat <- motifdat[motifid,]
  sapply(1:length(motifs),function(i) {
    motifs[[i]]$name <<- motifdat[i,'DBID.1']
    motifs[[i]]$tags <<- motifdat[i,]
  })
  motifs <- lapply(motifs,function(x) do.call(PWMatrix,x))
  motifs <- do.call(PWMatrixList,motifs)
  names(motifs) <- ID(motifs)
  return(motifs)
}

subMeta <- function(dev,sel){
  dev <- dev[sel,]
  metadata(dev) <- metadata(dev)[sel,]
  return(dev)
}

zfilt <- function(dev,zCutoff) subMeta(dev,apply(
    dev@assays$data$z,1,function(x) any(x>opts$zCutoff)
  ))

