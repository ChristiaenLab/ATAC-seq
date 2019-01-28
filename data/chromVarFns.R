chromVarPlots <- function(out,group='conditiontime',date=Sys.Date()) system2(
  'Rscript',paste(
    'chromVarPlots.R --dir',out,#paste0(out,'/',date),
    '--group',group
  )
)

# # chromVarPlots(opts$out,groups,Sys.Date())
# chromVarPlotDate <- function(
#   date=Sys.Date(),group='conditiontime',out='CISBP/Ciona_intestinalis_2017_11_01_3_41_pm/chromVAR/'
# ) sapply(
#   list.dirs(out,recursive=F),
#   function(x) chromVarPlots(paste0(x,'/',date),group))
# # chromVarPlotDate(groups,'2017-11-03')

chromVarKmers <- function(
  out,k=6,cisbp='Ciona_intestinalis_2017_11_01_3_41_pm'
) system2(
  'Rscript',paste(
    'chromVarKmers.R --dir',out,'--k',as.character(k),'--cisbp',cisbp
  )
)

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

motifHeatmap <- function(dev,de.motif=NULL,file,path=Sys.Date(),zCutoff=F,col=NULL,diff_acc=NULL,fdr=1,...){
  require(ComplexHeatmap)
  require(chromVAR)
  if(!is.null(de.motif)) dev <- subMeta(dev,de.motif)
  if(zCutoff) dev <- zfilt(dev,zCutoff)
  if(!is.null(diff_acc)){
    dev <- dev[intersect(names(dev),row.names(diff_acc)[diff_acc$p_value_adjusted<fdr])]
    diff_acc <- diff_acc[row.names(dev),]
    logp <- -log10(diff_acc$p_value_adjusted)
  }
  mat <- dev@assays$data$z
  tf.family <- data.frame(
    Family=unlist(elementMetadata(dev)$Family_Name),stringsAsFactors = F
  )
  if(!is.null(col)){
    sapply(names(col$Family), function(x){
      tf.family$Family[grep(x,tf.family$Family)] <<- x
    })
    # mat <- mat[tf.family$Family%in%names(col$Family),]
    
  }
  tf.family$Family <- forcats::fct_relevel(
    tf.family$Family,
    "bZIP","SAND","Sox","Nuclear receptor","Forkhead",
    "GATA","bHLH","Ets","Homeodomain","T-box"
  )
  colnames(mat) <- dev@colData$Name
  row.names(mat) <- make.names(elementMetadata(dev)$DBID.1,T)
  ht <- Heatmap(mat,row_names_gp = gpar(cex=.3),column_names_gp = gpar(cex=.8))
  
  tab <- cbind(
    tf.family, as.data.frame(ht@matrix)
  )#[row_order(ht)[[1]],]
  
  
  dir.eps(file,path,append.date = F)
  if(!is.null(col)){
    col$Family[Filter(function(x) !x%in%names(col$Family),levels(tf.family$Family))] <- 'white'
    ha <- rowAnnotation(tf.family,col=col,...)
    ht <- ht+ha
  }
  if(!is.null(diff_acc)){
    ha <- HeatmapAnnotation(
      points=anno_barplot(
        logp,which = 'row',axis = T
      ),
      which = 'row',annotation_width = unit(2,'cm')
    )
    ht <- ht+ha
    tab <- cbind(mat,diff_acc,logp)
  }
  draw(ht)
  dev.off()
  dir.tab(
    tab,file,path,append.date = F
  )
  
  htsplit <- Heatmap(mat,row_names_gp = gpar(cex=.3),column_names_gp = gpar(cex=.8),split=tf.family)
  dir.eps(paste0(file,'Split'),path,append.date = F)
  draw(htsplit)
  dev.off()
}

plotDenovo <- function(de_novos,motifs,counts,out,subdir='denovoMotifs',dist.cutoff=F){
  require(chromVAR)
  require(BSgenome.Cintestinalis.KH.KH2013)
  require(TFBSTools)
  require(motifmatchr)
  require(ggplot2)
  require(ggseqlogo)
  require(gridBase)
  require(gridExtra)
  # load(paste0(out,'/deviations.Rdata'))
  dist_to_known <- pwmDistance(de_novos, motifs)
  denovo_dist <- apply(dist_to_known$dist,1,min)
  which_match <- apply(dist_to_known$dist,1,which.min)
  closest_match <- do.call(PWMatrixList,motifs@listData[which_match])
  
  if(dist.cutoff){
    sel <- denovo_dist<dist.cutoff
    de_novos <- de_novos[sel]
    closest_match <- closest_match[sel]
    denovo_dist <- denovo_dist[sel]
  }
  
  denovo_matches <- matchMotifs(
    de_novos, counts, genome = BSgenome.Cintestinalis.KH.KH2013)
  
  # computing deviations
  denovo_dev <- computeDeviations(
    object = counts,  annotations = denovo_matches
  )
  
  denovo_cor <- mapply(
    cor,
    as.data.frame(t(denovo_dev@assays$data$deviations)),
    as.data.frame(t(dev@assays$data$deviations[which_match,]))
  )
  denovo_var <- computeVariability(denovo_dev)
  denovo_dat <- data.frame(
    name=denovo_var$name,
    denovo_dist,denovo_cor,
    denovo_var=denovo_var$variability,
    known_var=dev_var[which_match,"variability"]
  )
  
  out <- addSubdir(out,subdir)
  dir.tab(
    cbind(denovo_dat,do.call(rbind,lapply(tags(closest_match),unlist))),
    'closestMatch.txt',out,append.date = F
  )
  zones <- rbind(
    1:5+length(de_novos)*2+3,
    cbind(
      replicate(2,as.vector(rbind(
        1:length(de_novos),
        (length(de_novos)+1):(length(de_novos)*2)
      ))),
      length(de_novos)*2+1,
      length(de_novos)*2+2,
      length(de_novos)*2+3
    )
  )
  zones[,1] <- zones[,1]+length(de_novos)*2+3
  
  rm_axis <- theme(
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), 
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(), 
    axis.ticks.y=element_blank() 
  )
  rm_yaxis <- theme(
    axis.title.x=element_blank(), 
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(), 
    axis.ticks.y=element_blank() 
  )
  denovo_bar <- ggplot(denovo_dat)+rm_yaxis+coord_flip(expand=F)+scale_x_reverse()+geom_col()
  
  dist_bar <- denovo_bar+aes(x=as.numeric(name),y=denovo_dist)
  cor_bar <- denovo_bar+aes(x=as.numeric(name),y=denovo_cor)
  var_bar <- ggplot(
    data= reshape::melt(denovo_dat[,c("name","known_var","denovo_var")])
  )+aes(
    x=as.numeric(name),y=value,fill=variable
  )+geom_col(position = 'dodge' )+coord_flip(expand=F)+rm_yaxis+scale_x_reverse()
  
  logos <- lapply(Matrix(de_novos), function(x) ggseqlogo(x)+rm_axis)
  logos <- append(logos,lapply(Matrix(closest_match), function(x) ggseqlogo(x)+rm_axis))
  
  
  width <- c(2,4,4,4,6)
  height <- rep(1,nrow(zones))
  plot.denovo <- function(...) grid.arrange(
    ...,dist_bar,cor_bar,var_bar,
    layout_matrix=zones,
    heights=unit(height,"inches"),
    widths=unit(width,"inches"))
  
  dir.pdf('denovo',out,width=sum(width),height=sum(height),append.date = F)
  do.call(plot.denovo,logos)
  dev.off()
  
  lapply(names(de_novos),function(x) {
    dir.tab(
      t(Matrix(de_novos[[x]])),
      x,out,append.date = F
    )
    dir.pdf(x,out,width=sum(width),height=sum(height),append.date = F)
    seqLogo::seqLogo(Matrix(de_novos)[[x]])
    dev.off()
    ggsave(
      paste0(out,x,'Mismatch.eps'),
      plotKmerMismatch(de_novos[[x]]@tags$seed,kmer_cov),
      'eps',width=5,height=4
    )
  })
}