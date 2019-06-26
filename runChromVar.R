library(TFBSTools)
library(BiocParallel)
library(DBI)

# number of threads to use
register(MulticoreParam(20))

source('data/chromVarFns.R')
source('data/dirfns.R')
source("data/sqlfns.R")
source("data/getSelex.R")
source("data/getMotifs.R")

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ann <- getAnnotation(con)
peaksets <- getPeaksets(con)
scrna <- getScRNA(con)
# motifs <- getHomerMotifs("known.motifs")
motifs <- reduceMotifs(con,F,F,F)

mespPeaks <- ann$peaks[peaksets$mespDep]

expDesign <- dbReadTable(con,"ataclib",row.names="lib")
row.names(expDesign) <- paste0(sub('^X','',row.names(expDesign)),'_q30_rmdup_KhM0_sorted')
row.names(expDesign) <- paste0(row.names(expDesign),'.bam')

expDesign <- expDesign[
  expDesign$tissue=='B7.5'&expDesign$omit==0&expDesign$KO!=1,c("Name","condition","time")
]

mespDesign <- expDesign[
  expDesign$time%in%c('6hpf','10hpf'),
]

dev <- getChromVAR(
  expDesign,
  ann$peaks[Reduce(union,peaksets[c("timeDep","mespDep","handrDep")])],
  motifs
)

mespDev <- getChromVAR(mespDesign,mespPeaks,motifs)

denovoCardiacPeaks <- ann$peaks[
  unique(geneToPeak(con,scrna$denovoCardiac)$PeakID)
]
cardiacDev <- getChromVAR(expDesign,denovoCardiacPeaks,motifs)

denovoASMPeaks <- ann$peaks[
  unique(geneToPeak(con,scrna$denovoASM)$PeakID)
]
asmDev <- getChromVAR(expDesign,denovoASMPeaks,motifs)

asmCardiacDev <- getChromVAR(expDesign,ann$peaks[
  unique(geneToPeak(con,Reduce(
    union,scrna[c('denovoASM','denovoCardiac')]
  )))
])

save(
  dev,mespDev,cardiacDev,asmDev,asmCardiacDev,
  file = mkdate('chromVarOut','Rdata')
)

mapply(
  dbWriteTable,
  c('mapk10chromVAR','denovoCardiacChomVAR','denovoAsmChromVAR'),
  lapply(list(dev,mespDev,cardiacDev,asmDev),chromVarTable),
  MoreArgs = list(conn=con)
)