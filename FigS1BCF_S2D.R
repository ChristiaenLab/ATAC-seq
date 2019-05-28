# Figs. S1B,C,F,S2D
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(optparse)
library(DBI)
library(lattice)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
source('data/dirfns.R')
source("data/GRfns.R")
source("data/DESeqFns.R")
source('data/sqlfns.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
bulkGS <- getBulkRNA(con)
peakGeneAnnotation <- getAnnotation(con)
peaksets <- getPeaksets(con)
peakome <- peakGeneAnnotation$peaks

# Fig. S1
DEgenes <- intersect(union(unlist(scrna[-11]),unlist(bulkGS)),names(peakGeneAnnotation$geneToPeak))

de.daAnn <- mergeGenePeak(con,DEgenes,Reduce(union,peaksets))

de.daAnn <- list(
  geneToPeak=split(de.daAnn$PeakID,de.daAnn$GeneID),
  peakToGene=split(de.daAnn$GeneID,de.daAnn$PeakID)
)
DE.DAgeneLenBinom <- geneSizeBinom(
  names(de.daAnn$geneToPeak),de.daAnn,peakGeneAnnotation$genewindow[names(de.daAnn$geneToPeak)]
)
daAnn <- peakToGene(con,Reduce(union,peaksets))
daAnn <- list(
  geneToPeak=split(daAnn$PeakID,daAnn$GeneID),
  peakToGene=split(daAnn$GeneID,daAnn$PeakID)
)
DAgeneLenBinom <- geneSizeBinom(
  names(daAnn$geneToPeak),daAnn,peakGeneAnnotation$genewindow[names(daAnn$geneToPeak)]
)

featcount <- apply(peakGeneAnnotation$features[-1:-3],2,sum)

genomeSize <- sum(seqlengths(BSgenome.Cintestinalis.KH.JoinedScaffold))
genomeCoverage <- sum(width(peakome))/genomeSize

pfeat <- sapply(
  genomefeat, function(x) sum(width(x))
)/genomeSize

binom.feat <- mapply(
  binom.test,
  featcount,
  length(peakome),
  pfeat,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

peakCoverage <- sapply(genomefeat,function(x) reduce(
  pintersect(findOverlapPairs(peakome,reduce(unlist(x)))),
  ignore.strand=T
))

peakomebyfeat <- sapply(peakCoverage,function(x) sum(width(x)))/sum(width(peakome))
featbypeakome <- mapply(function(x,y) sum(width(x))/sum(width(y)),x=peakCoverage,y=genomefeat)


de.bulk.scrna <- Reduce(union,append(scrna,bulkGS))
de.bulk.scrna.peak <- Reduce(union,peakGeneAnnotation$geneToPeak[de.bulk.scrna])
de.count <- apply(peakGeneAnnotation$features[de.bulk.scrna.peak,-1:-3],2,sum)
de.binom.feat <- mapply(
  binom.test,
  de.count,
  length(de.bulk.scrna.peak),
  pfeat,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

tsc.count <- apply(peakGeneAnnotation$features[peakGeneAnnotation$features$TSC,-1:-3],2,sum)
tsc.binom.feat <- mapply(
  binom.test,
  tsc.count,
  featcount,
  sum(peakGeneAnnotation$features$TSC)/length(peakome),
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)


peakGeneFeat <- data.frame(
  npeaks=featcount,
  featInGenome=pfeat,
  featCoverageByPeakome=featbypeakome,
  peakomeCoverageByFeature=peakomebyfeat,
  peaksInFeature=sapply(binom.feat,'[[','estimate'),
  lowerConf=sapply(binom.feat,'[[','conf.int')[1,],
  upperConf=sapply(binom.feat,'[[','conf.int')[2,],
  p=sapply(binom.feat,'[[','p.value'),
  DEpeaks=de.count,
  DEpeaksInFeature=sapply(de.binom.feat,'[[','estimate'),
  DElowerConf=sapply(de.binom.feat,'[[','conf.int')[1,],
  DEupperConf=sapply(de.binom.feat,'[[','conf.int')[2,],
  DEp=sapply(de.binom.feat,'[[','p.value'),
  TSCpeaks=tsc.count,
  TSCpeaksInFeature=sapply(tsc.binom.feat,'[[','estimate'),
  TSClowerConf=sapply(tsc.binom.feat,'[[','conf.int')[1,],
  TSCupperConf=sapply(tsc.binom.feat,'[[','conf.int')[2,],
  TSCp=sapply(tsc.binom.feat,'[[','p.value')
)

# Fig. S1B
dir.pdf('percentpeakome')
barplot(
  t(peakGeneFeat[,c('featCoverageByPeakome','peakomeCoverageByFeature')])*100,
  las=3,beside = T,legend.text = c(
    "% feature covered by peakome","%peakome covered by feature"
  ),args.legend = list(x='topleft'),
  ylim = c(0,max(peakGeneFeat[-nrow(peakGeneFeat),c('featCoverageByPeakome','peakomeCoverageByFeature')]*100))
)
dev.off()

dir.tab(peakGeneFeat,'peakGeneCoverage')

# Fig. S1C
dir.eps('peakome.genome.feat')
barplot(
  t(peakGeneFeat[-nrow(peakGeneFeat),c(
    'DEpeaksInFeature','peaksInFeature','featInGenome'
  )])*100,beside = T, 
  # names.arg = c(
  #   'Promoter 1kb',"Promoter 0.5kb","TSS","5' UTR","CDS","Intron","3' UTR","TTS","Intergenic"
  # ),
  las=2,legend.text = c(
    '% DE peaks overlapping feature','% peaks overlapping feature','% genome'
  ),args.legend = list(x='topleft')
)
mapply(
  segments,
  seq(1.5,by=4,length.out = nrow(peakGeneFeat)),peakGeneFeat$DElowerConf*100,
  seq(1.5,by=4,length.out = nrow(peakGeneFeat)),peakGeneFeat$DEupperConf*100
)
mapply(
  segments,
  seq(2.5,by=4,length.out = nrow(peakGeneFeat)),peakGeneFeat$lowerConf*100,
  seq(2.5,by=4,length.out = nrow(peakGeneFeat)),peakGeneFeat$upperConf*100
)
dev.off()

# Fig. S1F
dir.eps('peaksize')
hist(width(peakGeneAnnotation$peaks),1000,xlim=c(0,1000),border = 'gray',col='gray',ylim=c(0,8000))
abline(v=c(mean(width(peakGeneAnnotation$peaks)),median(width(peakGeneAnnotation$peaks))),lty=c(1,2))
dev.off()

# Fig. S2D
dir.eps('peakomeTSCfeat')
barplot(
  t(peakGeneFeat[c(-nrow(peakGeneFeat)+1),'TSCpeaksInFeature',drop=F])*100,beside = T, 
  las=2,ylim = c(0,max(peakGeneFeat[c(-nrow(peakGeneFeat)+1),'TSCupperConf']*100))
)
mapply(
  segments,
  seq(1.5,by=2,length.out = nrow(peakGeneFeat)-1),
  peakGeneFeat[c(-nrow(peakGeneFeat)+1),'TSClowerConf']*100,
  seq(1.5,by=2,length.out = nrow(peakGeneFeat)-1),
  peakGeneFeat[c(-nrow(peakGeneFeat)+1),'TSCupperConf']*100
)
dev.off()
