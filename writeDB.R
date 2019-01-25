library(DBI)
library(GenomicRanges)
# library(reshape2)
source('data/sqlfns.R')
source('data/dirfns.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

# appends list of data.frames as in reshape2::melt
melt.rename <- function(dat,value='value',...){
  # require(reshape2)
  # res <- melt(dat)
  res <- sapply(dat,as.data.frame,simplify = F)
  res <- mapply(
    cbind,dat,names(dat),stringsAsFactors=F,SIMPLIFY = F
  )
  res <- lapply(res, function(x) setNames(x,c(names(x)[-length(x)],value)))
  res <- do.call(rbind,res)
  return(res)
}

dbWriteGenes <- function(...) dbWriteKey(...,foreign = 'gene_name(GeneID)')
dbWritePeaks <- function(...) dbWriteKey(...,foreign = 'peakfeature(PeakID)')

# this date will depend on when annotatePeaks.R was run
load('2018-12-05/peakome/peakGeneAnnotation.Rdata')

dbWriteKey(con,'peakfeature',peakGeneAnnotation$features,primary = 'PeakID',row.names = "PeakID")

genedat <- lrtab("dat/genes/",read.delim,'txt',quote='',stringsAsFactors=F)
genedat$ATM_genes_from_ANISEED <- genedat$ATM_genes_from_ANISEED[!duplicated(genedat$ATM_genes_from_ANISEED),]

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
# mapply(
#   dbWriteGenes,
#   name=names(scrna),
#   value=scrna,
#   MoreArgs = list(conn=con)
# )

genedat <- lrtab("dat/genes/",read.csv,'csv',stringsAsFactors=F)

names(genedat$mesenchyme20hpf)[1] <- "UniqueNAME"
# dbWriteTable(con,'mesenchyme20hpf',genedat$mesenchyme20hpf,overwrite=T)
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

gene.peak <- reshape2::melt(peakGeneAnnotation$peakToGene)
names(gene.peak) <- c("GeneID","PeakID")
dbWriteKey(con,'geneToPeak',gene.peak,foreign = c("gene_name(GeneID)","peakfeature(PeakID)"))

 # bulkRNAseq

foxf <- lrtab('dat/rnaseq/foxf/',read.csv,'\\.csv',stringsAsFactors=F)
# fold changes are reversed in these files
foxf$FoxF10hpf_LacZ10hpf$log2FoldChange <- -foxf$FoxF10hpf_LacZ10hpf$log2FoldChange
foxf$LacZ10hpf_Ngn10hpf$log2FoldChange <- -foxf$LacZ10hpf_Ngn10hpf$log2FoldChange

# foxf <- do.call(rbind,mapply(function(x,y) cbind(
#   as.data.frame(x),comparison=y,stringsAsFactors=F
# ),foxf,names(foxf),SIMPLIFY = F))
dbWriteGenes(con,'foxf_rnaseq',melt.rename(foxf,'comparison'))

handr <- lrtab('dat/rnaseq/handr/',read.csv,'\\.csv')
dbWriteGenes(con,'handr_rnaseq',melt.rename(handr,'comparison'))
# mapply(
#   dbWriteGenes,
#   name=names(rna),
#   value=rna,
#   MoreArgs = list(conn=con,row.names="GeneID")
# )

# ATACseq

peakcoord <- paste(
  seqnames(peakGeneAnnotation$peaks),
  start(peakGeneAnnotation$peaks)-1,
  end(peakGeneAnnotation$peaks),sep = ':'
)
peak.id.coord <- data.frame(
  PeakID=row.names(peakGeneAnnotation$features),
  row.names = peakcoord
)

metadat <- lrtab('dat/meta',read.delim,'txt',quote="'")
mapply(dbWriteTable,name=names(metadat),metadat,MoreArgs = list(conn=con,overwrite=T))

metadat <- lrtab('dat/meta/',read.csv,'csv',row.names=1)

mapply(
  dbWritePeaks,
  name=names(metadat),
  metadat,
  MoreArgs = list(conn=con)
)

peakdat <- lrtab('dat/peaks/',read.csv,'csv',row.names=1)
peakdat <- lapply(peakdat,function(x){
  row.names(x) <- peak.id.coord[row.names(x),]
  return(x)
})

mapply(
  dbWritePeaks,
  name=names(peakdat),
  peakdat,
  MoreArgs = list(conn=con,row.names="PeakID",overwrite=T)
)

ataclib <- dbReadTable(con,'ataclib')

counts <- lrtab(
  'counts.new/',
  read.table
)
dat <- do.call(cbind,sapply(counts,'[',7))
row.names(dat) <- peak.id.coord[do.call(function(...)paste(...,sep = ':'),counts[[1]][,1:3]),]
colnames(dat) <- make.names(sub('_counts_new.V7','',colnames(dat)))
dat <- as.data.frame(dat)
dat <- dat[,ataclib[,1]]
ataclib$reads <- apply(dat,2,sum)
dbWriteTable(con,'ataclib',ataclib,overwrite=T)

dbWritePeaks(con,'atacreads',dat,row.names = 'PeakID',overwrite=T)

atac <- Reduce(function(x,y){
  load(paste0(
    # the ATAC-seq directory corresponding to the date when runDESeq.R was run 
    '2018-10-02/peakomenocooks/', y, '/res.Rdata'
  ))
  x[[y]] <- res
  return(x)
  # the subdirectories corresponding to the DESeq2 models
}, c('foxf10','cardiac15','handr','lacz','gfp'), init = NULL)

names(atac$cardiac15) <- 'condition_handr_dnFGFR_15hpf_vs_control'
atac <- Reduce(append,atac)

atac <- lapply(atac,function(x){
  row.names(x) <- peak.id.coord[row.names(x),]
  return(x)
})

atac <- sapply(atac,function(x) cbind(PeakID=row.names(x),as.data.frame(x),stringsAsFactors=F),simplify = F)
dbWritePeaks(con,'atacseq',melt.rename(atac,'comparison'))

# atac <- do.call(rbind,mapply(function(x,y) cbind(
#   PeakID=row.names(x),as.data.frame(x),comparison=y,stringsAsFactors=F
# ),atac,names(atac),SIMPLIFY = F))

# mapply(
#   dbWritePeaks,
#   name=names(atac),
#   value=atac,
#   MoreArgs = list(conn=con,row.names='PeakID',overwrite=T)
# )

# dbWriteRownamesAs(con,'cisbp_matches',as.data.frame(as.matrix(motifMatches(peakmatches))),row.names = "PeakID")
# dbCreateTable(con,'cisbp_matches',colnames(motifMatches(peakmatches)))
# apply(
#   as.data.frame(as.matrix(motifMatches(peakmatches))),1,
#   function(x) dbWriteRownamesAs(con,'cisbp_matches',x,row.names = "PeakID",append=T)
# )

