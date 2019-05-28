source("data/DESeqFns.R")
promoterRna <- do.call(data.frame,lrtab('promoter_count',pattern = 'count',row.names=1))
names(promoterRna) <- sub('_promoter.*','',list.files('promoter_count/','count'))
# rnaDesign <- data.frame(
#   condition=c(
#     'FoxF_KO',rep("FoxF_camRas",7),"FoxF_KO",rep("FoxF_LacZ",8),rep("handr_dnFGFR",8),
#     rep("handr_LacZ",8),rep("LacZ_KO",2),rep("Ngn_KO",2)
#   ),
#   time=c(
#     "10hpf","12hpf",rep("15hpf",2),rep("18hpf",2),rep("20hpf",2),"10hpf",
#     rep(c(rep("12hpf",2),rep("15hpf",2),rep("18hpf",2),rep("20hpf",2)),3),
#     rep("10hpf",4)
#   ),
#   control=grepl('lacz',names(promoterRna)),
#   row.names = names(promoterRna)
# )
# rnaDesign$condtime <- paste0(rnaDesign$condition,'_',rnaDesign$time)
# write.csv(rnaDesign,'promoter_count/expDesign.csv')

de.promoterRna <- get.dds(promoterRna,data.frame(condition=c(
  rep('FoxF_KO_10hpf',2),rep('LacZ_10hpf',2),rep('Ngn_KO_10hpf',2)
),row.names = names(promoterRna)),~condition,list(
  c("condition","FoxF_KO_10hpf","LacZ_10hpf"),c("condition","FoxF_KO_10hpf","Ngn_KO_10hpf")
))

rnaDesign <- read.csv('promoter_count/expDesign.csv',row.names = 1)

de.promoterRna <- get.dds(
  promoterRna,
  rnaDesign[rnaDesign$experiment!='KO',],
  ~0+condition+experiment+time,
  'RNApromoterCount',
  reduced=~0+condition+experiment,
  test="LRT"
)

load('2018-07-03/RNApromoterCount/dds.Rdata')
de.promoterRna <- apply(
  rbind('time',combn(c('12hpf','15hpf','18hpf','20hpf'),2)),
  2,
  function(x) results(dds,x,format = 'DataFrame',alpha = .05,parallel = T)
)

promoterGtf <- import('peakome/promoters.gtf')
DEpromoter <- promoterGtf[promoterGtf$gene_id%in%row.names(de.promoterRna[[1]])[de.promoterRna[[1]]$padj<.05]]
DEpromoter <- split(DEpromoter,DEpromoter$gene_id)
DEpromoter <- unlist(reduce(DEpromoter))
DEpromoterPeaks <- findOverlaps(DEpromoter,peakGeneAnnotation$peaks)
findOverlaps(peakGeneAnnotation$peaks[to(DEpromoterPeaks)],tsc)

dir.export(DEpromoter,'DEpromoter')

library(motifmatchr)
library(GenomicRanges)
library(rtracklayer)
library(TFBSTools)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)
library(DBI)

source('data/sqlfns.R')
source('data/motifHyperFns.R')
source('data/chromVarFns.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

atac <- getAtac(con)
peakGeneAnnotation <- getAnnotation(con)

peakfeature <- peakGeneAnnotation$features
peaks <- peakGeneAnnotation$peaks
peakfeature <- peakfeature[,c(-1:-3,-13)]
sapply(1:length(peakfeature),function(x) peakfeature[,x] <<- as.logical(peakfeature[,x]))

motifs <- getCisbpMotifs('Ciona_intestinalis_2017_11_01_3_41_pm')
familyID <- sapply(tags(motifs),'[[',"Family_Name")
familyID <- sub('Homeodomain',"HD",familyID)
familyID <- sub('Nuclear receptor',"NR",familyID)
familyID <- sub("Forkhead","FH",familyID)
names(familyID) <- sapply(tags(motifs),'[[','TF_Name')
DEpromoterPeaks <- findOverlaps(DEpromoter,peaks,ignore.strand=T)

tscPeaks <- findOverlaps(tsc,peaks)

promoterTSSsel <- apply(peakfeature[,c('TSS','promoter500','promoter1k')],1,any)
tscsel <- intersect(which(promoterTSSsel),to(tscPeaks))
promotersel <- apply(peakfeature[,c('promoter500','promoter1k')],1,any)


peakmatches <- matchMotifs(motifs,peaks,genome = BSgenome.Cintestinalis.KH.JoinedScaffold,bg='genome')
featsel <- cbind(
  promoter=promotersel,peakfeature,
  # intergenic=!apply(peakGeneAnnotation$features,1,any),
  TSC=1:length(peakGeneAnnotation$peaks)%in%tscsel,
  DE=1:length(peakGeneAnnotation$peaks)%in%unique(to(DEpromoterPeaks))
)

promoterHyper <- motifHyper(peakmatches,featsel)

dir.tab(promoterHyper,'peakHyper')

corHyper <- data.frame(sapply(
  c(14,17:25),function(x) sapply(
    c(14,17:25), function(y) cor(
      promoterHyper[,x],promoterHyper[,y],method='spearman'
    )
  )
),row.names = colnames(promoterHyper)[c(14,17:25)])
names(corHyper) <- colnames(promoterHyper)[c(14,17:25)]

dir.eps('motifHyperCor')
Heatmap(corHyper)
dev.off()

dir.eps('haberleMotifHyper')
Heatmap(
  -log2(promoterHyper)[names(consensus)[c(-3,-4,-20)],c(26,29:37)],
  col = colorRamp2(c(0,-log2(0.05),20),c('blue','white','red')),
  width=unit(7,'cm'),
  heatmap_legend_param = list(title='log2(FDR)',color_bar="discrete")
)+rowAnnotation(
  consensus=row_anno_text(
    consensus[c(-3,-4,-20)],just='left',offset=unit(.1,"npc")
  ),
  annotation_width=unit(4,'cm')
)
dev.off()

dir.tab(t(matchHyper(sapply(atac,is.sig,p=.1),featsel)),'LFC1','atacFeatHyper')
dir.tab(t(matchHyper(sapply(atac,is.sig,lfc=.75,p=.1),featsel)), "LFC0.75","atacFeatHyper")
dir.tab(t(matchHyper(sapply(atac,is.sig,lfc=.5,p=.1),featsel)), "LFC0.5","atacFeatHyper")
dir.tab(t(matchHyper(sapply(atac,is.sig,lfc=0,p=.1),featsel)), "LFC0","atacFeatHyper")
# 
# dir.eps('barHyper',width=14)
# barplotHyper(cbind(
#   Any=res$condition_FoxF_KO_vs_control$padj<p.cutoff[2],
#   mapply(is.sig,res[c(
#     "condition_mesp_MekMut_vs_control","condition_mesp_dnFGFR_vs_control",
#     'condition_FoxF_KO_vs_control',
#     "condition_handr_dnFGFR_vs_control","condition_handr_MekMut_vs_control",
#     'time_6hpf_vs_10hpf','time_10hpf_18hpf'
#   )],lfc=c(.5,rep(.8,7)))
# ),featsel[,c(-8,-9,-11,-12)],xlim=c(0,130),col=colorRamp2(
#   c(1,round((ncol(featsel)-3)/2),ncol(featsel)-3),c('blue','white','red')
# )(1:(ncol(featsel)-3)))
# dev.off()
dir.eps('barHyper2',width=14)
barplotHyper(
  mapply(is.sig,atac[c(2:4,7:9,11,12)],lfc=.5,p=.05,tail='both'),
  featsel[,-9:-10],xlim=c(0,130),col=colorRamp2(
    c(1,round((ncol(featsel)-2)/2),ncol(featsel)-2),c('blue','white','red')
  )(1:(ncol(featsel)-2))
)
dev.off()

dir.eps('barHyper',width=14)
barplotHyper(
  sapply(
    list(
      TVC=tvcAcc,ATM=atmAcc,FoxF=closedFoxf,
      Cardiac=heartAcc,ASM=asmAcc,
      openAt6=open6,closedAt6=closed6,closedAt18=closed18,openAt18=open18
    ),
    function(x) names(peakGeneAnnotation$peaks)%in%x
  ),
  featsel[,-9:-10],xlim=c(0,130),col=colorRamp2(
    c(1,round((ncol(featsel)-1)/2),ncol(featsel)-1),c('blue','white','red')
  )(1:(ncol(featsel)-1))
)
dev.off()


venn <- list(
  promoterPeaks=names(peakGeneAnnotation$peaks)[promotersel],
  TSSseq=as.character(1:length(tsc)),
  DEpromoterPeaks=names(DEpromoter)
)
venn$TSSseq[tscPeaks@from] <- names(peakGeneAnnotation$peaks)[tscPeaks@to]
venn$DEpromoterPeaks[from(DEpromoterPeaks)] <- names(peakGeneAnnotation$peaks)[to(DEpromoterPeaks)]
de.tsc.promoter <- findOverlaps(tsc,DEpromoter)
venn$TSSseq[from(de.tsc.promoter)] <- venn$DEpromoterPeaks[to(de.tsc.promoter)]
VennDiagram::venn.diagram(venn,mkdate('TSSseq.venn','tiff'))

barplot(-log10(motifHyper$fdr.promoter))

promoterRna <- do.call(data.frame,lrtab('promoter_count',row.names=1,pattern = '.txt'))
countsRna <- do.call(data.frame,lrtab('~/ciona/star',pattern = 'count.txt',row.names=1))
sel <- intersect(row.names(promoterRna),row.names(countsRna))
promoterRna <- promoterRna[sel,]
countsRna <- countsRna[sel,]
names(promoterRna) <- sub('_promoter.*','',list.files('promoter_count/','.txt'))
names(countsRna) <- names(promoterRna)
mapply(cor,promoterRna,countsRna)

de.promoterRna <- get.dds(promoterRna,data.frame(condition=c(
  rep('FoxF_KO_10hpf',2),rep('LacZ_10hpf',2),rep('Ngn_KO_10hpf',2)
),row.names = names(promoterRna)),~condition,list(
  c("condition","FoxF_KO_10hpf","LacZ_10hpf"),c("condition","FoxF_KO_10hpf","Ngn_KO_10hpf")
))

ebfpromoter <- unique(peakGeneAnnotation$peaks[findOverlaps(import('ebfpromoter.bed'),peakGeneAnnotation$peaks)@to])
ebfmatch <- matchMotifs(
  motifs[familyID%in%c('Ets','T-box','ohler','EBF',"E-box")],
  unique(ebfpromoter),bg='subject',
  genome = BSgenome.Cintestinalis.KH.JoinedScaffold
)
ebfmotifs <- apply(
  motifMatches(ebfmatch), 1, function(x) cbind(
    motif=ebfmatch$name[x],
    Family=familyID[familyID%in%c('Ets','T-box','ohler','EBF',"E-box")][x]
  )
)
ebfmotifs <- mapply(cbind,peak=c('putative_promoter','10hpf-B','10hpf-A','pinky'),ebfmotifs)
ebfmotifs <- do.call(rbind,ebfmotifs)
dir.tab(ebfmotifs,'EBFpromoterMotifs')

peakmatches <- matchMotifs(motifs,peakGeneAnnotation$peaks,genome = BSgenome.Cintestinalis.KH.JoinedScaffold,bg='subject')


tscOverlap <- subsetByOverlaps(tsc,peakGeneAnnotation$peaks)

promoterct <- apply(motifMatches(peakmatches)[promotersel,],2,sum)
tscct <- apply(motifMatches(peakmatches)[tscsel,],2,sum)
DEpromoterct <- apply(motifMatches(peakmatches)[unique(to(DEpromoterPeaks)),],2,sum)

cdsct <- apply(motifMatches(peakmatches)[peakGeneAnnotation$features[,"CDS"],],2,sum)
intergenicct <- apply(motifMatches(peakmatches)[!apply(peakGeneAnnotation$features,1,any),],2,sum)
intronct <- apply(motifMatches(peakmatches)[peakGeneAnnotation$features[,"intron"],],2,sum)
ttsct <- apply(motifMatches(peakmatches)[peakGeneAnnotation$features[,"TTS"],],2,sum)
fputrct <- apply(motifMatches(peakmatches)[peakGeneAnnotation$features[,"five_prime_utr"],],2,sum)
tputrct <- apply(motifMatches(peakmatches)[peakGeneAnnotation$features[,"three_prime_utr"],],2,sum)

promoterHyper <- mapply(
  phyper,q=promoterct,m=peakct,k=sum(promotersel),n=n.peakct,lower.tail=F
)
tscHyper <- mapply(
  phyper,q=tscct,m=peakct,k=length(tscsel),n=n.peakct,lower.tail=F
)
deHyper <- mapply(
  phyper,q=DEpromoterct,m=peakct,
  k=length(unique(to(DEpromoterPeaks))),
  n=n.peakct,
  lower.tail=F
)


# tscOverlapMatch <- matchMotifs(
#   motifs, tscOverlap, genome=BSgenome.Cintestinalis.KH.JoinedScaffold, bg='genome'
# )
# tscOverlapCt <- apply(motifMatches(tscOverlapMatch),2,sum)
# tscOverlapHyper <- mapply(
#   phyper,q=tscOverlapCt,m=peakct,k=length(tscOverlap), n=n.peakct, lower.tail=F
# )

peakomeHyper <- data.frame(
  motif=peakmatches@colData$name,n.bg=peakct,
  n.promoter=promoterct,p.promoter=promoterHyper,fdr.promoter=p.adjust(promoterHyper,'fdr'),
  n.tsc=tscct,p.tsc=tscHyper,fdr.tsc=p.adjust(tscHyper,'fdr'),
  n.DE=DEpromoterct,p.DE=deHyper,fdr.DE=p.adjust(deHyper,'fdr')
  # n.tscOverlap=tscOverlapCt,p.tscOverlap=tscOverlapHyper,fdr.tscOverlap=p.adjust(tscOverlapHyper,'fdr')
)
dir.tab(peakomeHyper,'promoterHyper')
