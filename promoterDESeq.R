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
