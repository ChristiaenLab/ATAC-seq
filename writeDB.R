# initialize the database used for all subsequent analysis

library(DBI)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/sqlfns.R')
source('data/dirfns.R')
source('data/GRfns.R')

# initialize database
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

# read gene data
genedat <- lrtab("dat/genes/",read.delim,'txt',quote='',stringsAsFactors=F)
genedat$ATM_genes_from_ANISEED <- genedat$ATM_genes_from_ANISEED[!duplicated(genedat$ATM_genes_from_ANISEED),]
kh2013 <- getFeatures('KH.KHGene.2013.gff3',genedat$gene_name)

# TSS-seq
tsc <- read.csv('dat/peaks/tsc.csv')
tsc$GeneID <- paste0("KH2013:",tsc$GeneID)
dbWriteGenes(con,'tsc',tsc)

tsc <- GRanges(
  Rle(tsc$Chr),
  IRanges(tsc$Start,tsc$End),
  Rle(tsc$Strand),
  GeneID=tsc$GeneID,
  Rep.TSS=tsc$Rep.TSS,
  feature=tsc$Location
)

# annotate peaks
peakome <- import('peaks.subtracted.bed')
peakome <- peakome[width(peakome)>50]

promoterAnn <- findOverlaps(peakome,kh2013$promoterGene)

genomefeat <- append(
  kh2013$features,
  list(TSC=tsc,genome=GRanges(seqinfo(BSgenome.Cintestinalis.KH.KH2013)))
)
genomefeat <- do.call(GRangesList,lapply(genomefeat,function(x) reduce(unlist(x))))
genomefeat <- sapply(genomefeat,trim)

peakGeneAnnotation <- annotatePeaks(peakome,kh2013$genebody,features = genomefeat)
peakGeneAnnotation$promoterAnn <- promoterAnn


dbWriteKey(con,'peakfeature',peakGeneAnnotation$features,primary = 'PeakID',row.names = "PeakID")

gene.peak <- reshape2::melt(peakGeneAnnotation$peakToGene)
names(gene.peak) <- c("GeneID","PeakID")
dbWriteKey(con,'geneToPeak',gene.peak,foreign = c("gene_name(GeneID)","peakfeature(PeakID)"))


genebody <- data.frame(
  GeneID=names(peakGeneAnnotation$genes),
  chr=seqnames(peakGeneAnnotation$genes),
  start=start(peakGeneAnnotation$genes)+10000,
  end=end(peakGeneAnnotation$genes)-10000,
  stringsAsFactors = F
)
add.genes <- genebody$GeneID[!genebody$GeneID%in%genedat$gene_name$GeneID]
genedat$gene_name[add.genes,] <- rep(add.genes,2)
genedat$gene_name <- merge(genebody,genedat$gene_name)
dbWriteKey(con,"gene_name",genedat$gene_name,primary = "GeneID")

mapply(
  dbWriteGenes,
  name=names(genedat),
  genedat,
  MoreArgs = list(conn=con)
)

# scRNAseq
scrna <- lrtab('dat/scrna/',read.csv,stringsAsFactors=F)

scrna$TVCP <- as.data.frame(rbind(as.matrix(scrna$TVCP),cbind(
  c("KH2013:KH.L71.1","KH2013:KH.C10.370","KH2013:KH.C1.404","KH2013:KH.C1.279"),
  dbReadTable(con,'gene_name',row.names="GeneID")[c(
    "KH2013:KH.L71.1","KH2013:KH.C10.370","KH2013:KH.C1.404","KH2013:KH.C1.279"
  ),"UniqueNAME"],"TVC-specific","Primed"
)))

dbWriteGenes(con,'scrna',melt.rename(scrna,'geneset'))

ebfdat <- lrtab('dat/ebfdat/',read.csv,stringsAsFactors=F)
dbWriteGenes(con,'ebfdat',melt.rename(ebfdat,'geneset'))

genedat <- lrtab("dat/genes/",read.csv,'csv',stringsAsFactors=F)

names(genedat$mesenchyme20hpf)[1] <- "UniqueNAME"
dbWriteKey(con,'mesenchyme20hpf',genedat$mesenchyme20hpf,foreign = "gene_name(UniqueNAME)")

genedat$mesenchyme20hpf <- do.call(data.frame,append(
  dbGetQuery(
    con,
    'SELECT GeneID FROM mesenchyme20hpf LEFT JOIN gene_name ON mesenchyme20hpf.UniqueNAME=gene_name.UniqueNAME'
  ),
  genedat$mesenchyme20hpf
))

mapply(
  dbWriteGenes,
  name=names(genedat),
  genedat,
  MoreArgs = list(conn=con)
)

 # bulkRNAseq

foxf <- lrtab('dat/rnaseq/foxf/',read.csv,'\\.csv',stringsAsFactors=F)
# fold changes are reversed in these files
foxf$FoxF10hpf_LacZ10hpf$log2FoldChange <- -foxf$FoxF10hpf_LacZ10hpf$log2FoldChange
foxf$LacZ10hpf_Ngn10hpf$log2FoldChange <- -foxf$LacZ10hpf_Ngn10hpf$log2FoldChange

dbWriteGenes(con,'foxf_rnaseq',melt.rename(foxf,'comparison'))

handr <- lrtab('dat/rnaseq/handr/',read.csv,'\\.csv')
dbWriteGenes(con,'handr_rnaseq',melt.rename(handr,'comparison'))

# ATAC-seq metadata

metadat <- lrtab('dat/meta',read.delim,'txt',quote="'")
mapply(dbWriteTable,name=names(metadat),metadat,MoreArgs = list(conn=con,overwrite=T))

ataclib <- dbReadTable(con,'ataclib')

#generate PeakIDs
peakcoord <- paste(
  seqnames(peakGeneAnnotation$peaks),
  start(peakGeneAnnotation$peaks)-1,
  end(peakGeneAnnotation$peaks),sep = ':'
)
peak.id.coord <- data.frame(
  PeakID=row.names(peakGeneAnnotation$features),
  row.names = peakcoord
)

# ATAC-seq counts
counts <- lrtab(
  'counts.new/',
  read.table
)
dat <- do.call(cbind,sapply(counts,'[',7))
row.names(dat) <- peak.id.coord[do.call(function(...)paste(...,sep = ':'),counts[[1]][,1:3]),]
colnames(dat) <- make.names(sub('_counts_new.V7','',colnames(dat)))
dat <- as.data.frame(dat)
dat <- dat[,ataclib[,1]]
dbWritePeaks(con,'atacreads',dat,row.names = 'PeakID',overwrite=T)

# add reads per library to metadata
ataclib$reads <- apply(dat,2,sum)
dbWriteTable(con,'ataclib',ataclib,overwrite=T)


