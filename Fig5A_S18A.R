library(motifmatchr)
library(TFBSTools)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)

source('data/plotMotifs.R')
source('data/sqlfns.R')
source('data/chromVarFns.R')
source('data/motifHyperFns.R')
source('data/dirfns.R')
source('data/corHeatmap.R')
source('data/selex.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
ann <- getAnnotation(con)

known.motifs <- getHomerMotifs('known.motifs')
matches <- matchMotifs(known.motifs,ann$peaks,BSgenome.Cintestinalis.KH.JoinedScaffold,'subject','matches')
tf.family <- sapply(tags(known.motifs),'[[',"Family_Name")

cisbp.motifs <- getCisbpMotifs()
names(cisbp.motifs) <- make.names(sapply(tags(cisbp.motifs),'[[',"DBID.1"),T)
cisbp.matches <- matchMotifs(cisbp.motifs,ann$peaks,BSgenome.Cintestinalis.KH.JoinedScaffold,'subject','matches')
cisbp.family <- sapply(tags(cisbp.motifs),'[[',"Family_Name")

mekmut.dnfgfr.18 <- getAtacLib(con,c('condition_handr_MekMut_vs_control','condition_handr_dnFGFR_vs_control'))

peaks <- mapply(
  intersect,
  list(
      primedCardiacPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoCardiacPeaks=setdiff(peaksets$open18,union(peaksets$tvcAcc,peaksets$closed6)),
      primedAsmPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoAsmPeaks=setdiff(peaksets$open18,union(peaksets$tvcAcc,peaksets$closed6))
  ),
  list(
    # names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.25|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .25
    ],
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.25|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .25
    ],
    # names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.25|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.25
    ],
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.25|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.25
    ]
  )
)

tmp <- mekmut.dnfgfr.18$condition_handr_MekMut_vs_control[,1]-mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control[,1]
peaks <- mapply(
  intersect,
  list(
      primedCardiacPeaks=setdiff(peaksets$tvcAcc,Reduce(union,peaksets[c("open6","open18")])),
      denovoCardiacPeaks=setdiff(peaksets$open18,Reduce(union,peaksets[c("tvcAcc","closed6","open6")])),
      primedAsmPeaks=setdiff(peaksets$tvcAcc,Reduce(union,peaksets[c("open6","open18")])),
      denovoAsmPeaks=setdiff(peaksets$open18,Reduce(union,peaksets[c("tvcAcc","closed6","open6")]))
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

sapply(dynamics,length)


selex.hocomoco <- which(
  sapply(tags(cisbp.motifs),'[[',"Motif_Type")%in%c("SELEX","HocoMoco")&
    !grepl("^M[Afvw][0-9]+",names(cisbp.motifs))|
    names(cisbp.motifs)%in%c("EBF","K562_GATA2_UChicago")
)

# Fig. 5A
peakFamilyHyper(
  dynamics,
  'selex_hocomoco_TVC_18Family',
  cisbp.family[selex.hocomoco],
  motifMatches(cisbp.matches)[,selex.hocomoco],
  p=.050,fdr=T,logOR = 0.00,maskOR=T
)
colnames(cardiacDev) <- colData(cardiacDev)$Name
colnames(asmDev) <- colData(asmDev)$Name
colnames(asmCardiacDev) <- colData(asmCardiacDev)$Name
motifHmap(
  con,
  lapply(dynamics,lHyper,motifMatches(matches)),
  deviationScores(asmCardiacDev[,-6])[sel,],
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'denovoCardiac',or=1,padj=F
)

motifHmap(
  con,
  lapply(dynamics,lHyper,motifMatches(matches)),
  deviationScores(dev)[sel,],
  motifs,
  Reduce(union,append(
    bulkGS[c(-3:-4)],
    scrna[c(-8,-9,-11)]
  )),
  'fig5a',or=1,padj=F
)

# Fig. S18A
peakHyper(
  dynamics,
  'selex_hocomoco_TVC_18',
  motifMatches(cisbp.matches)[,selex.hocomoco],
  cisbp.family[selex.hocomoco],
  p=.050,fdr=F,logOR = 1.70,maskOR=T,breaks = c(0,4)
)

# Table S13
tmp <- do.call(rbind,peak.gene)
tmp$peakset <- sub('\\.[0-9]+','',row.names(tmp))
tmp$geneset <- c(
  rep('denovoCardiac',sum(sapply(peak.gene,nrow)[1:2])),
  rep('denovoASM',sum(sapply(peak.gene,nrow)[3:4]))
)

dir.tab(tmp,'primed_denovo_cardiac_asm',row.names=F)

motifs <- getSelex(con)
selex.family <- sapply(tags(motifs),'[[',"TF.family")

matches <- matchMotifs(motifs,ann$peaks,BSgenome.Cintestinalis.KH.JoinedScaffold,'subject','matches')
peakHyper(
  dynamics,
  'selex',
  motifMatches(matches),
  selex.family,
  p=0.100,fdr=F,logOR = 1,maskOR=T
)

up15downFgfr <- readLines('up15downFgfr.txt')
up15downMek <- readLines('up15downMek.txt')

peaks <- mapply(
  intersect,
  list(
      primedCardiacPeaks=setdiff(peaksets$tvcAcc,peaksets$open15),
      denovoCardiacPeaks=setdiff(peaksets$open15,peaksets$tvcAcc),
      primedAsmPeaks=setdiff(peaksets$tvcAcc,peaksets$open15),
      denovoAsmPeaks=setdiff(peaksets$open15,peaksets$tvcAcc)
  ),
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

peak.gene <- mapply(
  mergeGenePeak2,
  peaks=peaks,
  genes=list(up15downMek,up15downMek,up15downFgfr,up15downFgfr),
  MoreArgs = list(con=con),
  SIMPLIFY = F
)

dynamics <- lapply(peak.gene,function(x) unique(x$PeakID))
sapply(dynamics,length)

peakHyper(
  dynamics,
  'selex12_15',
  motifMatches(matches),
  selex.family,
  p=0.100,fdr=F,logOR = 1,maskOR=T
)

peakHyper(
  dynamics,
  'selex_hocomoco_12_15',
  motifMatches(cisbp.matches)[,selex.hocomoco],
  cisbp.family[selex.hocomoco],
  p=.100,fdr=F,logOR = 1.00,maskOR=T,breaks = c(0,4)
)

tmp <- do.call(rbind,peak.gene)
tmp$peakset <- sub('\\.[0-9]+','',row.names(tmp))
tmp$geneset <- c(
  rep('up15downMek',sum(sapply(peak.gene,nrow)[1:2])),
  rep('up15downFgfr',sum(sapply(peak.gene,nrow)[3:4]))
)

dir.tab(tmp,'primed_denovo_cardiac_asm_15',row.names=F)

peakfeat <- apply(ann$features,2,function(x) row.names(ann$features)[x])
feat15 <- unlist(lapply(
  list(heart15=up15downMek,asm15=up15downFgfr),
  function(x) lapply(
    peakfeat,
    function(y) mergeGenePeak(con,x,y)$PeakID
  )
),F)

peakHyper(
  feat15,
  'selex_hocomoco_12_15feat',
  motifMatches(cisbp.matches)[,selex.hocomoco],
  cisbp.family[selex.hocomoco],
  p=.100,fdr=F,logOR = 1.00,maskOR=T,breaks = c(0,4)
)

peakHyper(
  feat15,
  'selex_12_15feat',
  motifMatches(matches),
  selex.family,
  p=.100,fdr=F,logOR = 1.00,maskOR=T,breaks = c(0,4)
)

peakHyper(
  dynamics,
  'homer',
  motifMatches(matches),
  tf.family,
  p=.050,fdr=F,logOR = 1.70,maskOR=T#,breaks = c(0,4)
)

dir.tab(matchHyper(motifMatches(matches),sapply(peaksets,function(x)row.names(matches)%in%x)),'knownMotifHyper')
peakHyper(
  dynamics,
  'selex_dynamics',
  motifMatches(selex.matches8mer),
  selex.family,
  p=.100,fdr=F,logOR = 1.00,maskOR=T
)
