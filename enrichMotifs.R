library(TFBSTools)
library(motifmatchr)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)
library(chromVAR)
library(SummarizedExperiment)

source('data/chromVarFns.R')
source('data/sqlfns.R')
source('data/getMotifs.R')
source('data/alignMotifs.R')
source('data/motifHmap.R')
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ann <- getAnnotation(con)
peaksets <- getPeaksets(con)

bg <- letterFrequency(Views(Cintestinalis,ann$peaks),c("A","C","G","T"))
bg <- apply(bg,2,sum)
bg <- bg/sum(bg)

motifs <- reduceMotifs(con,F,F,F)
motif.dup <- setNames(motifs,name(motifs))
alignMotifs('nkx2_3.fa',motif.dup,bg=bg)
alignMotifs('smurf.fa',motif.dup,bg=bg)
alignMotifs('hand.fa',motif.dup,bg=bg)
alignMotifs('handfull.fa',motif.dup,bg=bg)

matches <- matchMotifs(motifs,ann$peaks,Cintestinalis,bg=bg)
sel <- split(names(motifs),ID(motifs))

fn <- function(x){
  if(length(x)>1){
    pwms <- combn(x,2)
    pearson <- apply(pwms,2,function(x) PWMSimilarity(
      motifs[[x[1]]],
      motifs[[x[2]]],
      "Pearson"
    ))
    return(x[!x%in%pwms[2,pearson>.90]])
  } else return(x)
}

motif.cor <- lapply(sel, fn)

dupct <- lapply(sel,function(x) apply(motifMatches(matches)[,x,drop=F],2,sum))
sel <- sapply(dupct,function(x) names(x)[which.max(x)])

hyper <- lapply(append(peaksets[c(5,6,1,2)],setNames(
  lapply(peaksets[c('open6','closed6')],intersect,peaksets$tvcAcc),
  c("earlyTVC","lateTVC")
),3),lHyper,motifMatches(matches))


load('2019-06-26/chromVarOut.Rdata')
dev <- dev[,c(-6,-12:-15)]
colnames(dev) <- colData(dev)$Name
diff_acc <- differentialDeviations(dev,'condtime')

mespDev <- mespDev[,-6]
colnames(mespDev) <- colData(mespDev)$Name
diff_acc <- differentialDeviations(mespDev,'condtime')

sel <- split(names(motifs),ID(motifs))
sel <- sapply(sel,function(x) x[which.min(diff_acc[x,1])])
sel <- intersect(unlist(sel),row.names(dev))

devscore <- deviationScores(mespDev)
devscore <- devscore[unlist(sel),]

motifHmap(
  con,
  hyper,
  devscore,
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifs',or=0
)

cardiacDev <- cardiacDev[,-6]
colnames(mespDev) <- colData(mespDev)$Name
diff_acc <- differentialDeviations(mespDev,'condtime')

sel <- split(names(motifs),ID(motifs))
sel <- sapply(sel,function(x) x[which.min(diff_acc[x,1])])
sel <- intersect(unlist(sel),row.names(dev))

devscore <- deviationScores(mespDev)
devscore <- devscore[unlist(sel),]

motifHmap(
  con,
  hyper,
  devscore,
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifs',or=0
)