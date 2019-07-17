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
dir.tab(cbind(
  ID=names(motifs),
  TF_Name=name(motifs),
  GeneID=ID(motifs),
  do.call(rbind,lapply(tags(motifs),'[',c("Family_Name","Motif_Type","DBID.1")))
),'motifdat',row.names=F,quote=T)

motif.dup <- setNames(motifs,name(motifs))
alignMotifs('nkx2_3.fa',motif.dup,bg=bg)
alignMotifs('smurf.fa',motif.dup,bg=bg)
alignMotifs('hand.fa',motif.dup,bg=bg)
alignMotifs('handfull.fa',motif.dup,bg=bg)

ebfscore <- matchMotifs(motifs,ann$peaks[
    paste0('KhL24.',as.character(c(37:34,31,28,27)))
],Cintestinalis,bg=bg,out='scores')

ebfscore <- sapply(motifByGene,function(x) apply(motifScores(ebfscore)[,x,drop=F],1,max))
colnames(ebfscore) <- sapply(motifByGene,'[',1)

plotPeakMatches(
  ebfscore,
  'ebfscore',
  sapply(tags(motifs),'[[',"Family_Name")[colnames(ebfscore)]
)

t12 <- matchMotifs(
  motifs,
  c(setNames(GRanges("KhC7",IRanges(1974923,1975287)),'t12'),ann$peaks[c(
    # "KhC7.909",
    "KhC7.914")]),
  Cintestinalis,bg=bg,out='scores'
)
tbxscore <- sapply(motifByGene,function(x) apply(motifScores(t12)[,x,drop=F],1,max))
colnames(tbxscore) <- sapply(motifByGene,'[',1)

plotPeakMatches(
  tbxscore,
  "tbxmotifs",
  sapply(tags(motifs),'[[',"Family_Name")[colnames(tbxscore)]
)

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
  'DEmotifs',or=0
)

peaks <- list(
  denovoCardiacPeaks=setdiff(
    peaksets$heartAcc,
    intersect(peaksets$tvcAcc,peaksets$open18)
  ),
  denovoAsmPeaks=setdiff(
    peaksets$asmAcc,
    intersect(peaksets$tvcAcc,peaksets$open18)
  )
)
peaks <- append(peaks,list(
  primedCardiacPeaks=setdiff(
    intersect(peaksets$tvcAcc,peaksets$open18),
    Reduce(union,peaks)
  ),
  primedAsmPeaks=setdiff(
    intersect(peaksets$tvcAcc,peaksets$open18),
    Reduce(union,peaks)
  )
))

dynamics <- mapply(
  function(x,y) unique(x[
    from(findOverlaps(ann$peaks[x],ann$genes[y]))
  ]),
  peaks[c(3,1,4,2)],
  scrna[c("denovoCardiac","denovoCardiac","denovoASM","denovoASM")]
)

peaks <- list(
    primedCardiacPeaks=setdiff(
      intersect(
        peaksets$tvcAcc,
        union(peaksets$open18,peaksets$heartAcc)
      ),peaksets$asmAcc
    ),
    denovoCardiacPeaks=setdiff(
      setdiff(
        union(
          peaksets$heartAcc,peaksets$open18
        ),peaksets$tvcAcc
      ),peaksets$asmAcc
    ),
    primedAsmPeaks=setdiff(
      intersect(
        peaksets$tvcAcc,union(
          peaksets$open18,peaksets$asmAcc
        )
      ),peaksets$heartAcc
    ),
    denovoAsmPeaks=setdiff(
      setdiff(
        union(
          peaksets$asmAcc,peaksets$open18
        ),peaksets$tvcAcc
      ),peaksets$heartAcc
    )
)

sapply(peaks,length)

peaks <- list(
  denovoCardiacPeaks=row.names(atac[[1]])[
    do.call('>',lapply(atac[c(
      "condition_handr_dnFGFR_vs_control",
      "condition_handr_MekMut_vs_control"
    )],'[',,1))
  ],
  denovoAsmPeaks=row.names(atac[[1]])[
    do.call('<',lapply(atac[c(
      "condition_handr_dnFGFR_vs_control",
      "condition_handr_MekMut_vs_control"
    )],'[',,1))
  ]
)

peaks <- lapply(peaks,intersect,peaksets$open18)
peaks <- mapply(union,peaks,peaksets[c('heartAcc','asmAcc')])
peaks <- mapply(setdiff,peaks,peaksets[c('asmAcc','heartAcc')])
peaks <- lapply(peaks,setdiff,Reduce(union,peaksets[c("tvcAcc",'closed18')]))

peak.gene <- mapply(
  mergeGenePeak2,
  peaks=peaks[c(3,1,4,2)],
  genes=scrna[c("denovoCardiac","denovoCardiac","denovoASM","denovoASM")],
  MoreArgs = list(con=con),
  SIMPLIFY = F
)

dynamics <- lapply(peak.gene,function(x) unique(x$PeakID))
sapply(dynamics,length)

asmCardiacDev <- asmCardiacDev[,-6]
colnames(asmCardiacDev) <- colData(asmCardiacDev)$Name
asmCardiacDiff <- differentialDeviations(asmCardiacDev[
  !apply(deviationScores(asmCardiacDev),1,anyNA),
],'condtime')
asmCardiacDiff <- sapply(motifByGene,function(x) x[which.min(asmCardiacDiff[x,1])])
asmCardiacDiff <- intersect(unlist(asmCardiacDiff),row.names(asmCardiacDev))

tmp <- mekmut.dnfgfr.18$condition_handr_MekMut_vs_control[,1]-mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control[,1]
peaks <- mapply(
  intersect,
  list(
      # primedCardiacPeaks=setdiff(peaksets$tvcAcc,Reduce(union,peaksets[c("open6","open18")])),
      # denovoCardiacPeaks=setdiff(peaksets$open18,Reduce(union,peaksets[c("tvcAcc","closed6","open6")])),
      # primedAsmPeaks=setdiff(peaksets$tvcAcc,Reduce(union,peaksets[c("open6","open18")])),
      # denovoAsmPeaks=setdiff(peaksets$open18,Reduce(union,peaksets[c("tvcAcc","closed6","open6")]))
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
      # mergeGenePeak
      #   con,
      #   Reduce(union,append(
      #     bulkGS[c(-3:-4)],
      #     scrna[c(-8,-9,-11)]
      #   )),
      #   Reduce(union,peaksets[c('mespDep','timeDep','handrDep')])
      # )$PeakID
    ,]
  ),
  deviationScores(asmCardiacDev[unlist(motif.cor),]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifsAsmCardiacExcl',.05,or=.5,F
)

peakFamilyHyper(
  dynamics,'denovoFamily',
  sapply(tags(motifs),'[[',"Family_Name"),
  motifMatches(matches),
  logOR = 0.00,fdr=F
)

cardiacDev <- cardiacDev[,-6]
colnames(cardiacDev) <- colData(cardiacDev)$Name
cardiacDiff <- differentialDeviations(cardiacDev[
  !apply(deviationScores(cardiacDev),1,anyNA),
],'condtime')
cardiacDiff <- sapply(motifByGene,function(x) x[which.min(cardiacDiff[x,1])])
cardiacDiff <- intersect(unlist(cardiacDiff),row.names(cardiacDev))

motifHmap(
  con,
  lapply(
    dynamics[1:2],lHyper,motifMatches(matches)
  ),
  deviationScores(cardiacDev[cardiacDiff,]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifsCardiac',.05,or=1,F
)

asmDev <- asmDev[,-6]
colnames(asmDev) <- colData(asmDev)$Name
asmDiff <- differentialDeviations(asmDev[
  !apply(deviationScores(asmDev),1,anyNA),
],'condtime')
asmDiff <- sapply(motifByGene,function(x) x[which.min(asmDiff[x,1])])
asmDiff <- intersect(unlist(asmDiff),row.names(asmDev))

motifHmap(
  con,
  lapply(
    dynamics[3:4],lHyper,motifMatches(matches)
  ),
  deviationScores(asmDev[asmDiff,]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifsAsm',.05,or=1,F
)


# devscore <- deviationScores(mespDev)
# devscore <- devscore[unlist(sel),]

motifHmap(
  con,
  hyper,
  deviationScores(mespDev[sel,]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifs2',or=0
)

motifHmap(
  con,
  hyper <- lapply(
    peaksets[c(5,6,1,2)],
    lHyper,
    motifMatches(matches)[Reduce(union,peaksets[c('timeDep','mespDep','handrDep')]),]
  ),
  deviationScores(mespDev[sel,]),
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'DEmotifsVsDA',or=0
)

mespDiff <- sapply(motifByGene,function(x) x[which.min(mespDiff[x,1])])
mespDiff <- intersect(unlist(sel),row.names(dev))
