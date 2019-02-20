library(motifmatchr)
library(TFBSTools)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/chromVarFns.R')
source('data/dirfns.R')

cisbp.motifs <- getCisbpMotifs()
names(cisbp.motifs) <- make.names(sapply(tags(cisbp.motifs),'[[',"DBID.1"),T)
cisbp.family <- sapply(tags(cisbp.motifs),'[[',"Family_Name")

peaks <- c("KhL24:31630-32784","KhL24:33162-33770","KhL24:32318-32848","KhL24:33162-33779")
ranges <- GRanges(peaks)

cisbp.matches <- matchMotifs(
  cisbp.motifs[cisbp.family=="promoter"],
  ranges,
  BSgenome.Cintestinalis.KH.KH2013,'subject','positions'
)

nmotifs <- sapply(cisbp.matches,countOverlaps,query=ranges)
row.names(nmotifs) <- peaks

dir.tab(t(nmotifs),'ebf_promoter_motifs')
nmotifs <- nmotifs/width(ranges)*1000

dir.tab(t(nmotifs),'ebf_promoter_motifs_per_kb')
