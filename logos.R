# logos shown in Fig. 2A

library(ggplot2)
library(ggseqlogo)

source('data/chromVarFns.R')
source('data/dirfns.R')

homer.motifs <- getHomerMotifs('known.motifs')
sapply(
  homer.motifs[sapply(tags(homer.motifs),'[[',"DBID.1")%in%c(
    'MyoD',"MyoG","c-Jun-CRE","Foxf1","Nkx6.1","Smad4","RARg","Pax7","Tbx5","Gata4","Gata6","Otx2","CRX"
  )],
  function(x) dir.gg(
    ggseqlogo(Matrix(x))+ggtitle(tags(x)$DBID.1),
    tags(x)$DBID.1,ext='eps'
  )
)
