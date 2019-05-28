getSelex <- function(){
  require(TFBSTools)
  
  selex <- read.csv('Selex_seq_Cirobu_All/SELEX_project_TableS1_S3_Nitta_et_al_2019.csv',stringsAsFactors = F)
  selexGene <- read.csv("Selex_seq_Cirobu_All/SELEXgeneID.csv",stringsAsFactors = F)
  selex <- merge(selex,selexGene,by.x=1,by.y=1)
  selex <- selex[selex$Selex.Motif.Found=="Yes",c(
    'Primary.gene.name','Best.transcript.model',"Family.Classification","Best.Pfm.6mer","Best.Pfm.8mer"
  )]
  selex$Best.transcript.model <- paste0("KH2013:",sub('.v.*','',selex$Best.transcript.model))
  selex <- selex[!duplicated(selex[,-1]),]
  selex[duplicated(selex$Primary.gene.name),"Primary.gene.name"] <- paste0(
    selex[duplicated(selex$Primary.gene.name),"Primary.gene.name"],'v2'
  )
  
  profileMatrix8mer <- apply(
    selex[,c("Best.Pfm.8mer","Best.Pfm.6mer")],
    1,
    function(x) {
      if(x[1]==''){
        file <- x[2]
      }else{
        file <- x[1]
      }
      return(as.matrix(
        read.table(paste0("Selex_seq_Cirobu_All/pwms/",file),row.names=c("A","C","G","T"))
      ))
    }
  )
  
  profileMatrix8mer <- lapply(profileMatrix8mer,function(x) t(t(x)/apply(x,2,sum)))
  selex.pwm8mer <- mapply(
    PWMatrix,
    selex$Primary.gene.name,
    selex$Best.transcript.model,
    profileMatrix=profileMatrix8mer,
    tags=lapply(selex$Family.Classification,function(x) list(family=x))
  )
  selex.pwm8mer <- do.call(PWMatrixList,selex.pwm8mer)
  return(selex.pwm8mer)
} 