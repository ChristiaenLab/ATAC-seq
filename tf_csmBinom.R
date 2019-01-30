# Figs. 
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(optparse)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)
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

tfs <- dbReadTable(con,'TFs_from_GHOST',row.names="GeneID")
csm <- dbReadTable(con,'signaling_molecules_from_GHOST',row.names="GeneID")

DEgenes <- intersect(union(unlist(scrna[-11]),unlist(bulkGS)),names(peakGeneAnnotation$geneToPeak))

tf.csmBinom <- geneSizeBinom(
  list(TF=tfs,CSM=csm),peakGeneAnnotation,peakGeneAnnotation$genewindow
)
DEtf.csmBinom <- geneSizeBinom(
  list(TF=tfs,CSM=csm),peakGeneAnnotation,peakGeneAnnotation$genewindow[DEgenes]
)

# Fig. S3G
dir.eps('DEtf.csmBar')
barplot(
  t(DEtf.csmBinom[,c('null.value','estimate')])*100,
  beside = T, 
  # names.arg = c(
  #   'Promoter 1kb',"Promoter 0.5kb","TSS","5' UTR","CDS","Intron","3' UTR","TTS","Intergenic"
  # ),
  las=2,legend.text = c(
    '% DE gene coverage','% peaks overlapping DE genes'
  ),args.legend = list(x='topright'),
  ylim=c(0,max(DEtf.csmBinom[,"upperConf"])*100)
)
segments(
  seq(2.5,by=3,length.out = nrow(DEtf.csmBinom)),DEtf.csmBinom[,"lowerConf"]*100,
  seq(2.5,by=3,length.out = nrow(DEtf.csmBinom)),DEtf.csmBinom[,"upperConf"]*100
)
dev.off()

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


windowSize <- sum(width(reduce(peakGeneAnnotation$genewindow)))
tfSize <- sum(width(reduce(peakGeneAnnotation$genewindow[tfs])))
csmSize <- sum(width(reduce(peakGeneAnnotation$genewindow[csm])))

DEtfs <- intersect(tfs,DEgenes)
DEcsm <- intersect(csm,DEgenes)
DEwindowSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEgenes])))
DEtfSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEtfs])))
DEcsmSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEcsm])))

tf.csmBinom <- mapply(
  binom.test,
  c(
    length(unique(unlist(peakGeneAnnotation$geneToPeak[tfs]))),
    length(unique(unlist(peakGeneAnnotation$geneToPeak[csm]))),
    length(unique(unlist(peakGeneAnnotation$geneToPeak[DEgenes])))
  ),
  length(peakGeneAnnotation$peakToGene),
 c(tfSize,csmSize,DEwindowSize)/windowSize,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

genelenBinom <- mapply(
  binom.test,
  sapply(peakGeneAnnotation$geneToPeak,length),
  length(peakGeneAnnotation$peakToGene),
  width(peakGeneAnnotation$genewindow[names(peakGeneAnnotation$geneToPeak)])/windowSize,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

binomDat <- cbind(
  n=sapply(genelenBinom,'[[','statistic'),
  sapply(genelenBinom,'[[','p.value'),
  p.adjust(sapply(genelenBinom,'[[','p.value')),
  sapply(genelenBinom,'[[','null.value'),
  sapply(genelenBinom,'[[','estimate')
)


DEtf.csmBinom <- mapply(
  binom.test,
  c(
    length(unique(unlist(peakGeneAnnotation$geneToPeak[DEtfs]))),
    length(unique(unlist(peakGeneAnnotation$geneToPeak[DEcsm])))
  ),
  length(peakGeneAnnotation$peakToGene),
  c(DEtfSize,DEcsmSize)/DEwindowSize,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

DEgenelenBinom <- mapply(
  binom.test,
  sapply(peakGeneAnnotation$geneToPeak[DEgenes],length),
  length(peakGeneAnnotation$peakToGene),
  width(peakGeneAnnotation$genewindow[DEgenes])/windowSize,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)

DEbinomDat <- cbind(
  n=sapply(DEgenelenBinom,'[[','statistic'),
  p=sapply(DEgenelenBinom,'[[','p.value'),
  padj=p.adjust(sapply(DEgenelenBinom,'[[','p.value')),
  null.value=sapply(DEgenelenBinom,'[[','null.value'),
  estimate=sapply(DEgenelenBinom,'[[','estimate'),
  lowerConf=sapply(DEgenelenBinom,'[[','conf.int')[1,],
  upperConf=sapply(DEgenelenBinom,'[[','conf.int')[2,]
)
dir.tab(DEbinomDat,"DEgeneLenBinom")

featcount <- apply(peakGeneAnnotation$features[-1:-3],2,sum)

genomeSize <- sum(seqlengths(BSgenome.Cintestinalis.KH.KH2013))
genomeCoverage <- sum(width(peakome))/genomeSize

pfeat <- sapply(
  genomefeat, function(x) sum(width(x))
)/genomeSize#sum(width(reduce(unlist(genomefeat))))

binom.feat <- mapply(
  binom.test,
  featcount,
  length(peakome),
  pfeat,
  alternative='two.sided',conf.level=.99,SIMPLIFY = F
)
# mapply(pbinom,feathits,sum(feathits),pfeat,lower.tail=F,log.p=F)

peakCoverage <- sapply(genomefeat,function(x) reduce(
  pintersect(findOverlapPairs(peakome,reduce(unlist(x)))),
  ignore.strand=T
))

peakomebyfeat <- sapply(peakCoverage,function(x) sum(width(x)))/sum(width(peakome))
featbypeakome <- mapply(function(x,y) sum(width(x))/sum(width(y)),x=peakCoverage,y=genomefeat)
# featbypeakome['genome'] <- genomeCoverage


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
tmp <- split(peakGeneAnnotation$peaks,seqnames(peakGeneAnnotation$peaks))
strand(tmp) <- "+"

peakGC <- letterFrequency(Views(BSgenome.Cintestinalis.KH.KH2013,peakGeneAnnotation$peaks),"GC")

featGC <- apply(peakGeneAnnotation$features,2,function(x) sum(peakGC[x])/sum(width(peakGeneAnnotation$peaks[x])))

peakGC/width(peakGeneAnnotation$peaks)

peakomeStr <- extractTranscriptSeqs(BSgenome.Cintestinalis.KH.KH2013,tmp)
peakomeFreq <- apply(letterFrequency(peakomeStr,letters = c("A","C","T","G")),2,sum)
peakomeFreq <- peakomeFreq/sum(peakomeFreq)
peakomeGC <- sum(peakomeFreq[c("C","G")])
genomeFreq <- sapply(
  seqnames(BSgenome.Cintestinalis.KH.KH2013),
  function(x)letterFrequency(BSgenome.Cintestinalis.KH.KH2013[[x]],c("A","C","T","G"))
)
genomeFreq <- apply(genomeFreq,1,sum)
genomeFreq <- genomeFreq/sum(genomeFreq)
genomeGC <- sum(genomeFreq[c("C","G")])

peakomeGC <- getGC(peakGeneAnnotation$peaks)
featGC <- apply(peakGeneAnnotation$features,2,function(x) getGC(peakGeneAnnotation$peaks[x]))
featGC$intergenic <- getGC(peakGeneAnnotation$peaks[!apply(peakGeneAnnotation$features,1,any)])
cdsGC <- getGC(reduce(unlist(cds)))

ncGC <- getGC(gaps(reduce(unlist(GRangesList(reduce(unlist(cds)),peakGeneAnnotation$peaks)))))
intersect(gaps(reduce(unlist(cds))),peakGeneAnnotation$peaks)
ncpeakGC <- getGC(setdiff(peakGeneAnnotation$peaks,reduce(unlist(cds)),ignore.strand=T))

barplot(c(CDS=cdsGC,NC=ncGC,peakome=peakomeGC)*100,ylim = c(30,50))
dir.eps('GCcontent')
barchart(c(CDS=cdsGC,NC=ncGC,peakome=peakomeGC)*100,xlim=c(30,50),col='blue')
dev.off()

#Fig. S1D
peakGC <- letterFrequency(Views(BSgenome.Cintestinalis.KH.KH2013,peakGeneAnnotation$peaks),"GC")

featGC <- apply(peakGeneAnnotation$features[,-1:-3],2,function(x) sum(peakGC[x])/sum(width(peakGeneAnnotation$peaks[x])))

genomeFeatGC <- sapply(genomefeat,function(x) sum(
  letterFrequency(Views(BSgenome.Cintestinalis.KH.KH2013,x),"GC")
)/sum(width(x)))

dir.eps('GCfeat')
barplot(rbind(featGC,genomeFeatGC),las=2,beside = T)
dev.off()

# Fig. S1E
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

