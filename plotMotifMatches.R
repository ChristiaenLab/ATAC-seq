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
  bg='subject',
  out='scores'
)

# select only motifs in desired families
sel <- homer.family%in%c(
  "AP2EREBP","bHLH","bZIP","E2F","EBF","ETS","ETS:IRF","Forkhead,NR","Forkhead,bHLH","Forkhead","Homeobox","ILF","IRF","MAD","MADS","NR","NRF","POU,Homeobox,IR1","SacCer-Promoters",'Stat',"T-box","Zf"
)

# plot heatmap for logical matrix of matches
plotPeakMatches(
  # extract matrix from homer.matches
  motifScores(homer.matches)[
    # subset by desired peaks
    paste0('KhL24.',as.character(c(37:34,31,28,27))),
    # subset by families 
    sel
  ],
  # file name
  'ebfmotifFamilySelex',
  # split heatmap by families
  # this should also be subset by sel, otherwise groups will be wrong
  selex.family[sel],
  # don't reorder columns
  cluster_columns=F
)

t12 <- matchMotifs(
  known.motifs,
  setNames(GRanges("KhC7",IRanges(1974923,1975287)),'t12'),
  BSgenome.Cintestinalis.KH.JoinedScaffold,
  "subject",
  out='scores'
)
plotPeakMatches(
  rbind(
    t12=motifScores(t12),
    KhC7.914=motifScores(homer.matches)[c(
      # "KhC7.909",
      "KhC7.914"
    ),]
  )[,sel],
  "tbxmotifs",
  homer.family[sel],
  cluster_columns=F
)

plotPeakMatches(
  motifScores(homer.matches)[
    c("KhC5.1641","KhC4.137","KhC4.144"),
    sel
  ],
  "mmp21_lrp4_8",
  homer.family[sel],
  cluster_columns=F
)

nkx <- import('nkx2_3.fa','fasta')
names(nkx) <- c('Crobusta','Csavignyi')
nkx.motifs <- matchMotifs(
  known.motifs,
  nkx,
  # GRanges('KhC8',IRanges(4058582,4060915)),
  # Cintestinalis,
  out='positions')
cr.nkx <- unlist(do.call(IRangesList,lapply(nkx.motifs,'[[',1)))
cr.nkx <- GRanges('Crobusta',do.call(IRanges,as.data.frame(cr.nkx)),mcols(cr.nkx)$strand,score=mcols(cr.nkx)$score)
cr.nkx$seq <- as.character(Views(nkx[[1]],ranges(cr.nkx)))
cs.nkx <- unlist(do.call(IRangesList,lapply(nkx.motifs,'[[',2)))
cs.nkx <- GRanges('Csavignyi',do.call(IRanges,as.data.frame(cs.nkx)),mcols(cs.nkx)$strand,score=mcols(cs.nkx)$score)
cs.nkx$seq <- as.character(Views(nkx[[2]],ranges(cs.nkx)))
nkx.matches <- c(cr.nkx,cs.nkx)
nkx.matches$family <- homer.family[names(nkx.matches)]
export(nkx.matches,'nkxMatches.gff3')

plotPeakMatches(
  motifScores(matchMotifs(
    known.motifs,nkx,out='scores'
  ))[,sel],
  'nkxConserved',
  homer.family[sel],
  column_labels=c('Crobusta','Csavignyi')
)

smurf <- import('smurf.fa','fasta')
smurf.motifs <- matchMotifs(
  known.motifs,
  smurf,
  out='positions')
cr.smurf <- unlist(do.call(IRangesList,lapply(smurf.motifs,'[[',1)))
cr.smurf <- GRanges('Crobusta',do.call(IRanges,as.data.frame(cr.smurf)),mcols(cr.smurf)$strand,score=mcols(cr.smurf)$score)
cr.smurf$seq <- as.character(Views(smurf[[1]],ranges(cr.smurf)))
cs.smurf <- unlist(do.call(IRangesList,lapply(smurf.motifs,'[[',2)))
cs.smurf <- GRanges('Csavignyi',do.call(IRanges,as.data.frame(cs.smurf)),mcols(cs.smurf)$strand,score=mcols(cs.smurf)$score)
cs.smurf$seq <- as.character(Views(smurf[[2]],ranges(cs.smurf)))
smurf.matches <- c(cr.smurf,cs.smurf)
smurf.matches$family <- homer.family[names(smurf.matches)]
export(smurf.matches,'smurfMatches.gff3')

plotPeakMatches(
  motifScores(matchMotifs(
    known.motifs,smurf,out='scores'
  ))[,sel],
  'smurfConserved',
  homer.family[sel],
  column_labels=c('Crobusta','Csavignyi')
)

nkx.motifs <- matchMotifs(selex.pwm8mer,GRanges('KhC8',IRanges(4058582,4060915)),Cintestinalis,out='positions')
nkx.motifs <- unlist(nkx.motifs)
nkx.motifs$seq <- as.character(Views(Cintestinalis,nkx.motifs))
nkx.motifs$family <- selex.family[names(nkx.motifs)]
export(nkx.motifs,'nkxMotifsSelex.gff3')

hand.motifs <- matchMotifs(known.motifs,GRanges('KhC14',IRanges(1826145,1827438)),Cintestinalis,out='positions')
hand.motifs <- unlist(hand.motifs)
hand.motifs$seq <- as.character(Views(Cintestinalis,hand.motifs))
hand.motifs$family <- homer.family[names(hand.motifs)]
export(hand.motifs,'handMotifs.gff3')

smurf.motifs <- matchMotifs(known.motifs,GRanges('KhC4',IRanges(1259557,1262258)),Cintestinalis,out='positions')
smurf.motifs <- unlist(smurf.motifs)
smurf.motifs$seq <- as.character(Views(Cintestinalis,smurf.motifs))
smurf.motifs$family <- homer.family[names(smurf.motifs)]
export(smurf.motifs,'smurfMotifs.gff3')

distanceToNearest(GRanges(
  c("KhC1","KhC6","KhC4","KhC4","KhC3","KhC2","KhC5"),
  IRanges(
    c(3167894,1856686,1906285,1902488,106453,6992204,4070723),
    c(3167913,1856710,1906304,1902509,106472,6992223,4070742)
  )
),ann$peaks)
