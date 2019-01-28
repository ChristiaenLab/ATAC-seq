# usage: Rscript chromVarDeviations.R --peaks <bam file> --conditiontime "<levels in atac/DEseq/expDesign.csv$conditiontime>" --out <out file name> --cores <number of parallel threads>
library(optparse)
library(chromVAR)
library(BSgenome.Cintestinalis.KH.KH2013)
library(TFBSTools)
library(BiocParallel)
library(motifmatchr)
library(DBI)
source('chromVarFns.R')
source('data/dirfns.R')
source("data/sqlfns.R")
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')


expDesign <- dbReadTable(con,"ataclib",row.names="lib")

row.names(expDesign) <- paste0(sub('^X','',row.names(expDesign)),'_q30_rmdup_KhM0_sorted')
# row.names(expDesign)[c(1:3,12,13)] <- paste0(row.names(expDesign)[c(1:3,12,13)],'_merged')
row.names(expDesign) <- paste0(row.names(expDesign),'.bam')
# row.names(expDesign)[!sapply(row.names(expDesign),file.exists)] <- sub('sorted','sorted_merged',row.names(expDesign)[!sapply(row.names(expDesign),file.exists)])

opts <- list(
  make_option('--peaks',type='character',default = 'peaks.subtracted.bed'),
  make_option('--conditiontime',type='character',default = 'all'),#unique(expDesign$conditiontime)),
  make_option('--condition',type='character',default = unique(expDesign$condition)),
  make_option('--time',type='character',default = unique(expDesign$time)),
  make_option('--tissue',type='character',default = unique(expDesign$tissue)),
  make_option('--noKO','store_true',default = F),
  make_option('--out',type='character',default = 'peakome'),
  make_option('--cisbp',type='character',default = 'Ciona_intestinalis_2017_11_01_3_41_pm'),
  make_option('--sub','store_true',default = F),
  make_option('--resize','store_true',default = T),
  make_option('--cores',default = 20,type = 'integer')
)
opts <- parse_args(OptionParser(option_list = opts))
out <- paste('chromVAR','homer',opts$out,sep='/')

register(MulticoreParam(opts$cores))

if(opts$noKO) expDesign <- expDesign[!expDesign$KO,]
expDesign <- expDesign[
  expDesign$condition%in%unlist(strsplit(opts$condition,','))&
  expDesign$time%in%unlist(strsplit(opts$time,','))&
  expDesign$tissue%in%unlist(strsplit(opts$tissue,','))&
  !expDesign$omit,
  c('Name','condition','time','tissue')
]


expDesign$condtime <- apply(expDesign[,-1],1,paste,collapse='')
ncondtime <- table(expDesign$condtime)
design <- expDesign[expDesign$condtime%in%names(ncondtime)[ncondtime>1],]

#motifs <- getCisbpMotifs(opts$cisbp)
motifs <- getHomerMotifs("known.motifs")


peaks <- import(opts$peaks)#getPeaks(opts$peaks,sort_peaks = T)
seqlengths(peaks) <- seqlengths(BSgenome.Cintestinalis.KH.KH2013)[seqlevels(peaks)]

if(opts$resize){
peaks <- resize(peaks,width = 200,fix = 'center')
peaks <- trim(peaks)
peaks <- resize(peaks,width = 200,fix = 'end')
peaks <- trim(peaks)
peaks <- resize(peaks,width = 200,fix = 'start')
}

counts <- getCounts(
  paste0('/scratch/kaw504/atacCiona/atac/bam/',row.names(design)),
  peaks, paired =  TRUE,  format = "bam", colData = DataFrame(design))
counts <- addGCBias(
  counts, genome = BSgenome.Cintestinalis.KH.KH2013)

matches <- matchMotifs(
  motifs, counts, genome = BSgenome.Cintestinalis.KH.KH2013
)

# computing deviations
dev <- computeDeviations(
  object = counts,  annotations = matches)
elementMetadata(dev) <- cbind(elementMetadata(dev),as.data.frame(do.call(rbind, tags(motifs))))
dev_var <- computeVariability(dev)
#dev_DA <- differentialDeviations(dev)

save(dev,dev_var,counts,matches,motifs,file=mkdate('deviations','Rdata',out))
