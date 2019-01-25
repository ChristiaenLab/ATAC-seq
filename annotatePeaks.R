library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(optparse)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/dirfns.R')
source("data/GRfns.R")
source('data/sqlfns.R')

opts <- list(
  make_option('--peakome',type = 'character',default='peaks.subtracted.bed'),
  make_option('--out',type = 'character',default = 'peakome'),
  make_option('--gff',type = 'character',default = 'KH.KHGene.2013.gff3')
)
opts <- OptionParser(
  "usage: Rscript --peakome <peakome.bed> --out <outdir> --gff <genome.gff3>",
  opts
)
opts <- parse_args(opts)

peakome <- import(opts$peakome)
peakome <- peakome[width(peakome)>50]
path <- opts$out

kh2013 <- getFeatures(import(opts$gff))

promoterAnn <- findOverlaps(peakome,kh2013$promoterGene)

tsc <- read.csv('tsc.csv')
tsc$GeneID <- paste0("KH2013:",tsc$GeneID)
tsc <- GRanges(
  Rle(tsc$Chromosome.Scaffold),
  IRanges(tsc$Start,tsc$End),
  Rle(tsc$Strand),
  GeneID=tsc$GeneID,
  Rep.TSS=tsc$Rep.TSS,
  feature=tsc$Location
)
export(tsc,mkdate('tsc','bed'))

genomefeat <- append(
  kh2013$features,
  list(TSC=tsc,genome=GRanges(seqinfo(BSgenome.Cintestinalis.KH.KH2013)))
)
genomefeat <- do.call(GRangesList,lapply(genomefeat,function(x) reduce(unlist(x))))
genomefeat <- sapply(genomefeat,trim)

peakGeneAnnotation <- annotatePeaks(peakome,kh2013$genebody,features = genomefeat)
peakGeneAnnotation$promoterAnn <- promoterAnn
save(peakGeneAnnotation,file = mkdate('peakGeneAnnotation','Rdata',path))

