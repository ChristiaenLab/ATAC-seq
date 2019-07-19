library(TFBSTools)
library(motifmatchr)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)
library(SummarizedExperiment)

source('data/dirfns.R')
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
# Fig S10C
alignMotifs('nkx2_3.fa',motif.dup,bg=bg)
# Fig S12C
alignMotifs('hand.fa',motif.dup,bg=bg)
alignMotifs('handfull.fa',motif.dup,bg=bg)
# Fig S13C
alignMotifs('smurf.fa',motif.dup,bg=bg)

# Fig S19A
ebfscore <- matchMotifs(motifs,ann$peaks[
    paste0('KhL24.',as.character(c(37:34,31,28,27)))
],Cintestinalis,bg=bg,out='scores')

ebfscore <- sapply(motifByGene,function(x) apply(motifScores(ebfscore)[,x,drop=F],1,max))
colnames(ebfscore) <- sapply(motifByGene,'[',1)

plotPeakMatches(
  ebfscore,
  'FigS19A',
  sapply(tags(motifs),'[[',"Family_Name")[colnames(ebfscore)]
)

# Fig S20C
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
  "FigS20C",
  sapply(tags(motifs),'[[',"Family_Name")[colnames(tbxscore)]
)