library(TFBSTools)
library(motifmatchr)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)
library(chromVAR)
library(SummarizedExperiment)

source('data/dirfns.R')
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

# Fig 3AB
hyper <- lapply(append(peaksets[c(5,6,1,2)],setNames(
  lapply(peaksets[c('open6','closed6')],intersect,peaksets$tvcAcc),
  c("earlyTVC","lateTVC")
),3),lHyper,motifMatches(matches))

hyper <- lapply(peaksets[c(5,6,1,2)],lHyper,motifMatches(matches))
motifByGene <- split(names(motifs),ID(motifs))

load('2019-06-26/chromVarOut.Rdata')
dev <- dev[,-6]
colnames(dev) <- colData(dev)$Name
diff_acc <- differentialDeviations(dev,'condtime')
diff_acc <- sapply(motifByGene,function(x) x[which.min(diff_acc[x,1])])
diff_acc <- intersect(unlist(diff_acc),row.names(dev))

mespDev <- mespDev[,-6]
colnames(mespDev) <- colData(mespDev)$Name
mespDiff <- differentialDeviations(mespDev,'condtime')
mespDiff <- sapply(motifByGene,function(x) x[which.min(mespDiff[x,1])])
mespDiff <- intersect(unlist(mespDiff),row.names(mespDev))


devscore <- deviationScores(mespDev)
devscore <- devscore[unlist(sel),]

motifHmap(
  con,
  hyper,
  deviationScores(mespDev[mespDiff,]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'Fig3AB',or=0
)

tmp <- mekmut.dnfgfr.18$condition_handr_MekMut_vs_control[,1]-mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control[,1]
peaks <- mapply(
  intersect,
  list(
      primedCardiacPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoCardiacPeaks=setdiff(peaksets$open18,peaksets$tvcAcc),
      primedAsmPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoAsmPeaks=setdiff(peaksets$open18,peaksets$tvcAcc)
  ),
  list(
    row.names(mekmut.dnfgfr.18[[1]])[tmp<0],
    row.names(mekmut.dnfgfr.18[[1]])[tmp<0],
    row.names(mekmut.dnfgfr.18[[1]])[tmp>0],
    row.names(mekmut.dnfgfr.18[[1]])[tmp>0]
  )
)
sapply(peaks,length)

peak.gene <- mapply(
  mergeGenePeak2,
  peaks=peaks,
  genes=scrna[c("denovoCardiac","denovoCardiac","denovoASM","denovoASM")],
  MoreArgs = list(con=con),
  SIMPLIFY = F
)

dynamics <- lapply(peak.gene,function(x) unique(x$PeakID))
dynamics <- sapply(dynamics,setdiff,geneToPeak(con,Reduce(union,append(
  scrna[c("TVCP","STVC")],
  bulkGS[c("MAPK10activated","downreg6hpf")]
)))$PeakID)

dynamics <- mapply(
  setdiff,
  dynamics,
  list(
    geneToPeak(con, Reduce(union,append(scrna["ASM"],bulkGS["MAPK18activated"])))$PeakID,
    geneToPeak(con, Reduce(union,append(scrna["ASM"],bulkGS["MAPK18activated"])))$PeakID,
    geneToPeak(con, Reduce(union,append(scrna["Cardiac"],bulkGS["MAPK18inhibited"])))$PeakiD,
    geneToPeak(con, Reduce(union,append(scrna["Cardiac"],bulkGS["MAPK18inhibited"])))$PeakiD
  )
)

motifHmap(
  con,
  lapply(
    dynamics,lHyper,motifMatches(matches)[
    ,]
  ),
  deviationScores(asmCardiacDev[unlist(motif.cor),]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'FigS18A',.05,or=.5,F
)

