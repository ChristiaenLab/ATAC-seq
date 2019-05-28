library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(motifmatchr)
library(rtracklayer)
library(TFBSTools)
library(DBI)
library(GenomicRanges)

source('data/chromVarFns.R')
source("data/dirfns.R")
source("data/plotMotifs.R")
source("data/corHeatmap.R")
# read peaks from database
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
peakGeneAnnotation <- getAnnotation(con)

# read homer motifs from file
known.motifs <- getHomerMotifs('known.motifs')
# extract families from motifs
homer.family <- sapply(tags(known.motifs),'[[',"Family_Name")

# get matches for each peak
homer.matches <- matchMotifs(
  known.motifs,
  peakGeneAnnotation$peaks,
  BSgenome.Cintestinalis.KH.JoinedScaffold,
  bg='subject'
)

# select only motifs in desired families
sel <- homer.family%in%c(
  "AP2EREBP","bHLH","bZIP","E2F","EBF","ETS","ETS:IRF","Forkhead,NR","Forkhead,bHLH","Forkhead","Homeobox","ILF","IRF","MAD","MADS","NR","NRF","POU,Homeobox,IR1","SacCer-Promoters",'Stat',"T-box","Zf"
)

# plot heatmap for logical matrix of matches
plotPeakMatches(
  # extract matrix from homer.matches
  motifMatches(selex.matches8mer)[
    # subset by desired peaks
    paste0('KhL24.',as.character(c(37:34,31,28,27))),
    # subset by families 
    # sel
  ],
  # file name
  'ebfmotifFamily',
  # split heatmap by families
  # this should also be subset by sel, otherwise groups will be wrong
  selex.family,#[sel],
  # don't reorder columns
  cluster_columns=F
)

t12 <- matchMotifs(
  selex.pwm8mer,
  setNames(GRanges("KhC7",IRanges(1974923,1975287)),'t12'),
  BSgenome.Cintestinalis.KH.JoinedScaffold,
  "subject"
)
plotPeakMatches(
  rbind(
    t12=motifMatches(t12),
    KhC7.914=motifMatches(selex.matches8mer)[c(
      # "KhC7.909",
      "KhC7.914"
    ),]
  ),
  "tbxmotifs",
  selex.family,
  cluster_columns=F
)

plotPeakMatches(
  motifMatches(selex.matches8mer)[c(
    "KhC5.1641","KhC4.137","KhC4.144"
  ),],
  "mmp21_lrp4_8",
  selex.family,
  cluster_columns=F
)

