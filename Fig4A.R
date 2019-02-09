library(motifmatchr)
library(TFBSTools)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/plotMotifs.R')
source('data/sqlfns.R')
source('data/chromVarFns.R')
source('data/motifHyperFns.R')
source('data/dirfns.R')
source('data/corHeatmap.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
ann <- getAnnotation(con)

known.motifs <- getHomerMotifs('known.motifs')
matches <- matchMotifs(known.motifs,ann$peaks,BSgenome.Cintestinalis.KH.KH2013,'subject','matches')
tf.family <- sapply(tags(known.motifs),'[[',"Family_Name")

cisbp.motifs <- getCisbpMotifs()
names(cisbp.motifs) <- make.names(sapply(tags(cisbp.motifs),'[[',"DBID.1"),T)
cisbp.matches <- matchMotifs(cisbp.motifs,ann$peaks,BSgenome.Cintestinalis.KH.KH2013,'subject','matches')
cisbp.family <- sapply(tags(cisbp.motifs),'[[',"Family_Name")

mekmut.dnfgfr.18 <- getAtacLib(con,c('condition_handr_MekMut_vs_control','condition_handr_dnFGFR_vs_control'))

peaks <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoPeaks=setdiff(
        peaksets$open18,
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)
peaks <- mapply(
  intersect,
  peaks,
  list(
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.00|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .00
    ],
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.00|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.00
    ]
  )
)

sapply(peaks,length)

selex.hocomoco <- which(
  sapply(tags(cisbp.motifs),'[[',"Motif_Type")%in%c("SELEX","HocoMoco")&
    !grepl("^M[Afvw][0-9]+",names(cisbp.motifs))|
    names(cisbp.motifs)%in%c("EBF","K562_GATA2_UChicago")
)

peakFamilyHyper(
  peaks,
  'selex_hocomoco_TVC_18Family',
  cisbp.family[selex.hocomoco],
  motifMatches(cisbp.matches)[,selex.hocomoco],
  p=.050,fdr=T,logOR = 0.00,maskOR=T
)

peakHyper(
  peaks,
  'selex_hocomoco_TVC_18',
  motifMatches(cisbp.matches)[,selex.hocomoco],
  cisbp.family[selex.hocomoco],
  p=.050,fdr=F,logOR = 1.70,maskOR=T,breaks = c(0,4)
)

tmp <- stack(sapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(
        Reduce(union,peaksets[c("open6",'closed6','closed18')]),
        peaksets$open18
      ),
      denovoPeaks=peaksets$open18
    ),
    function(y) mergeGenePeak2(con,x,y)
  )
))
dir.tab(tmp[,-3],'primed_denovo_cardiac_asm',col.names=F,row.names=F)

