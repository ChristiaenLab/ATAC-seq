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

alignMotifs <- function(fasta,motifs,species=c("Crobusta","Csavignyi"),suffix=''){
  require(rtracklayer)
  require(motifmatchr)
  require(TFBSTools)
  dat <- import(fasta,'fasta')
  names(dat) <- species
  family <- unlist(lapply(tags(motifs),'[[',"Family_Name"))
  matches <- matchMotifs(
    motifs,
    dat,
    out='positions'
  )
  matchpos <- sapply(1:length(species),function(x){
    cr <- lapply(matches,'[[',x)
    sel <- sapply(cr,length)>0
    family <- family[sel]
    cr <- cr[sel]
    cr <- mapply(function(x,y){
      mcols(x)$family <- y
      return(x)
    },cr,family)
    cr <- do.call(IRangesList,cr)
    cr <- unlist(cr)
    cr <- GRanges(
      species[x],
      do.call(IRanges,as.data.frame(cr)),
      mcols(cr)$strand,
      score=mcols(cr)$score,
      family=mcols(cr)$family
    )
    cr$seq <- as.character(Views(dat[[x]],ranges(cr)))
    return(cr)
  })
  names(matchpos) <- species
  # cr <- unlist(do.call(IRangesList,lapply(matches,'[[',1)))
  # cr <- GRanges(species[1],do.call(IRanges,as.data.frame(cr)),mcols(cr)$strand,score=mcols(cr)$score)
  # cr$seq <- as.character(Views(dat[[1]],ranges(cr)))
  # cs <- unlist(do.call(IRangesList,lapply(matches,'[[',2)))
  # cs <- GRanges(species[2],do.call(IRanges,as.data.frame(cs)),mcols(cs)$strand,score=mcols(cs)$score)
  # cs$seq <- as.character(Views(dat[[2]],ranges(cs)))
  matches <- Reduce(c,matchpos)
  # matches$family <- family[names(matches)]
  file <- paste0(sub('\\..*$','',fasta),suffix)
  dir.export(matches,paste0(file,"Matches"),format = "gff3")
  
  x <- sapply(matchpos,function(x) sapply(
    split(x,names(x)),
    function(y) mcols(y)[which.max(y$score),c('score','family')]
  ))
  x <- lapply(x,do.call,what=rbind)
  
  x <- Reduce(function(i,j) merge(i,j,c("row.names",'family'),all=T),x)
  y <- setNames(x[,-1:-2],species)
  y <- as.matrix(y)
  row.names(y) <- x[,1]#sub('(.{,20}).*','\\1',x[,1])
  # y[is.na(y)] <- 0
  sel <- !apply(y,1,anyNA)
  # x <- motifScores(matchMotifs(
  #   motifs,dat,out='scores'
  # ))
  # row.names(x) <- species
  plotPeakMatches(
    t(y[sel,]),
    paste0(file,'Conserved'),
    x$family[sel]
  )
}

tmp <- cisbp.motifs
names(tmp) <- sapply(
  split(
    sub("KH2013:",'',sub('.*_','',cisbpDat$GeneName)),
    cisbpDat$Motif_ID
  ),
  paste,collapse=';'
)[ID(cisbp.motifs)]
tmp <- tmp[!is.na(names(tmp))]

alignMotifs('nkx2_3.fa',tmp,suffix = 'CISBP')
alignMotifs('smurf.fa',tmp,suffix = 'CISBP')
alignMotifs('hand.fa',tmp,suffix = 'CISBP')
alignMotifs('handfull.fa',tmp,suffix = 'CISBP')

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
