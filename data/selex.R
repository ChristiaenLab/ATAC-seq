getSelex <- function(con,x='Selex_seq_Cirobu_bestround'){
  selex <- read.csv(paste0(x,'/SELEX_project_TableS1_S3_Nitta_et_al_2019.csv'),stringsAsFactors = F)
  pwms <- lrtab(paste0(x,'/pfm_best round'),pattern = '1-6-1')
  pwms <- lapply(pwms,function(x) t(t(unname(x))/apply(x,2,sum)))

  names(pwms) <- sub('D.*','',sub('_.*','',names(pwms)))
  
  motifs <- selex[selex$Selex.clone.AC.code%in%names(pwms),1:11]
  motifs <- do.call(PWMatrixList,apply(
    motifs,1,
    function(x) PWMatrix(
      ID=x["Selex.clone.AC.code"],
      tags = as.list(x),
      profileMatrix = as.matrix(setRowNames(
        pwms[[x["Selex.clone.AC.code"]]],
        c("A","C","G","T")
      ))
    )
  ))
  genes <- paste0("KH2013:",sub('\\.v.*','',sapply(tags(motifs),'[[',"Best.transcript.model")))
  
  names(motifs) <- make.names(sub('.*_','',getGeneNames(con)[genes,]),T)
  return(motifs)
}