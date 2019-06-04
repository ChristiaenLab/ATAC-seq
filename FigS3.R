# Fig. S3B-G
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

DEgenes <- intersect(union(unlist(scrna[-11]),unlist(bulkGS)),names(peakGeneAnnotation$geneToPeak))

# Fig. S3B
dir.eps('genesPerPeak')
hist(sapply(peakGeneAnnotation$peakToGene,length),border = 'gray',col='gray')
dev.off()
# Fig. S3C
dir.eps('peaksPerGene')
hist(sapply(peakGeneAnnotation$geneToPeak,length),100,border = 'gray',col='gray')
dev.off()

tfs <- intersect(
  dbReadTable(con,'TFs_from_GHOST')$GeneID,
  row.names(peakGeneAnnotation$gene.names)
)
csm <- intersect(
  dbReadTable(con,'signaling_molecules_from_GHOST')$GeneID,
  row.names(peakGeneAnnotation$gene.names)
)

windowSize <- sum(width(reduce(peakGeneAnnotation$genewindow)))
tfSize <- sum(width(reduce(peakGeneAnnotation$genewindow[tfs])))
csmSize <- sum(width(reduce(peakGeneAnnotation$genewindow[csm])))

DEtfs <- intersect(tfs,DEgenes)
DEcsm <- intersect(csm,DEgenes)
DEwindowSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEgenes])))
DEtfSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEtfs])))
DEcsmSize <- sum(width(reduce(peakGeneAnnotation$genewindow[DEcsm])))
# expressed <- Reduce(union,append(sapply(scrna,row.names),bulkGS))

peak.size <- data.frame(
  peaklen=width(peakGeneAnnotation$peaks),
  type='all',row.names = names(peakGeneAnnotation$peaks),stringsAsFactors=F
)
peak.size[unique(unlist(peakGeneAnnotation$geneToPeak[tfs])),'type'] <- "TF"
peak.size[unique(unlist(peakGeneAnnotation$geneToPeak[csm])),'type'] <- "CSM"

gene.size.npeak <- data.frame(
  bp = width(peakGeneAnnotation$genes[names(peakGeneAnnotation$geneToPeak)]),
  npeak = sapply(peakGeneAnnotation$geneToPeak,length),
  type='any',
  row.names = names(peakGeneAnnotation$geneToPeak),
  stringsAsFactors = F
)
gene.size.npeak$log2bp <- log2(gene.size.npeak$bp)
gene.size.npeak[tfs,'type'] <- "TF"
gene.size.npeak[csm,'type'] <- "CSM"
gene.size.npeak$peakPerKb <- gene.size.npeak$npeak/gene.size.npeak$bp*1000

# Fig. S3D
dir.eps('log2geneSize.npeak')
plot(gene.size.npeak[,c('log2bp','npeak')],pch=19,col='lightgray',cex=.8)
points(gene.size.npeak[tfs,c('log2bp','npeak')],pch=19,col='red',cex=.8)
points(gene.size.npeak[csm,c('log2bp','npeak')],pch=19,col='blue',cex=.8)
dev.off()

# Fig. S3E
dir.eps('DElog2geneSize.npeak')
plot(gene.size.npeak[,c('log2bp','npeak')],pch=19,col='lightgray',cex=.8)
points(gene.size.npeak[DEtfs,c('log2bp','npeak')],pch=19,col='red',cex=.8)
points(gene.size.npeak[DEcsm,c('log2bp','npeak')],pch=19,col='blue',cex=.8)
dev.off()

# Fig S3F
dir.eps('peakPerKb')
densityplot(~peakPerKb,gene.size.npeak,groups=type,auto.key = T,xlim = c(0,20),pch=NA)
dev.off()

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

# Fig. S3G
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

tf.csmBinom <- geneSizeBinom(
  list(TF=tfs,CSM=csm),peakGeneAnnotation,peakGeneAnnotation$genewindow
)
DEtf.csmBinom <- geneSizeBinom(
  list(TF=tfs,CSM=csm),peakGeneAnnotation,peakGeneAnnotation$genewindow[DEgenes]
)

dir.eps('DEtf.csmBar')
barplot(
  t(DEtf.csmBinom[,c('null.value','estimate')])*100,
  beside = T, 
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
