library(motifmatchr)
library(TFBSTools)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/chromVarFns.R')
source('data/dirfns.R')
source('data/plotMotifs.R')
source('data/corHeatmap.R')

cisbp.motifs <- getCisbpMotifs()
names(cisbp.motifs) <- make.names(sapply(tags(cisbp.motifs),'[[',"DBID.1"),T)
cisbp.family <- sapply(tags(cisbp.motifs),'[[',"Family_Name")

peaks <- c("KhL24:31630-32784","KhL24:33162-33770","KhL24:32318-32848","KhL24:33162-33779")
ranges <- GRanges(peaks)

cisbp.matches <- matchMotifs(
  cisbp.motifs,
  setNames(ranges,peaks),
  BSgenome.Cintestinalis.KH.KH2013,'subject','positions'
)

nmotifs <- sapply(cisbp.matches,countOverlaps,query=ranges)
row.names(nmotifs) <- peaks

dir.tab(t(nmotifs),'ebf_promoter_motifs')
nmotifs <- nmotifs/width(ranges)*1000

dir.tab(t(nmotifs),'ebf_promoter_motifs_per_kb')

known.motifs <- getHomerMotifs('known.motifs')
family <- sapply(tags(known.motifs),'[[',"Family_Name")

matches <- matchMotifs(
  known.motifs,
  setNames(ranges,c('a','b','c','d')),
  BSgenome.Cintestinalis.KH.KH2013,
  'subject','matches'
)

plotPeakMatches(motifMatches(matches),'ebfmatches',family)

library(ggseqlogo)
ggseqlogo(Matrix(motifs$FOXF1.2))


cisbp.matches <- cisbp.matches[sapply(cisbp.matches,length)>0]

tmp <- lapply(cisbp.matches,Views,subject=BSgenome.Cintestinalis.KH.KH2013)

ranges <- sapply(tmp,granges)
ranges <- do.call(GRangesList,ranges)
ranges <- unlist(ranges)

mcols(ranges)$seq <- unlist(sapply(tmp,as.character))

library(rtracklayer)
dir.export(ranges,'cisbpMatch',path='..',format='gff3')
