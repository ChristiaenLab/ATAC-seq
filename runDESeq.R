source("data/DESeqFns.R")
source('data/dirfns.R')

expDesign <- read.delim('dat/meta/ataclib.txt',row.names = 1)
expDesign <- subset(expDesign,!omit)

# read counts from htseq
counts <- lrtab(
  'counts/',
  read.table
)

dat <- do.call(cbind,sapply(counts,'[',4))
row.names(dat) <- do.call(function(...)paste(...,sep = ':'),counts[[1]][,1:3])
colnames(dat) <- make.names(sub('_counts_new.V4','',colnames(dat)))
dat <- as.data.frame(dat)
dat <- dat[,row.names(expDesign)]
# dir.tab(dat,'counts','peakome')

# get peak lengths
len <- do.call('-',counts[[1]][,3:2])
# remove peaks shorter than 50bp
filt <- dat[len>50,]
# remove replicates with fewer than 500000 reads
filt <- filt[,apply(filt,2,sum)>500000]
# reorder expDesign by column order of filt
expDesign <- expDesign[names(filt),]

# removes the bottom 5% of peaks, not used
rm.quant <- function(dat,q=.05) dat[rowMeans(dat)>quantile(rowMeans(dat),q),]

# tissue-dependent model
design.gfp <- subset(expDesign,condition=='control'&tissue!='Ef1a'&time!='6hpf'&!KO)
dat.gfp <- filt[,row.names(design.gfp)]
get.dds(
  filt,design.gfp,~0+time*tissue+method,
  path='peakomenocooks/gfp',
  test="LRT",reduced=~0+time+tissue+method,cooksCutoff=F,alpha=.05

)

# subset peaks by B7.5 lineage
b75expDesign <- subset(expDesign,tissue=='B7.5')
b75 <- filt[,row.names(b75expDesign)]
# b75 <- rm.quant(b75)
filt <- filt[row.names(b75),]

# 10hpf model
design10 <- subset(b75expDesign,time=='10hpf')
get.dds(
  filt[,row.names(design10)],design10,~0+condition+KO,
  path='peakomenocooks/cond10',
  comparisons=list(
   c('condition','FoxF_KO',"control"),
    c('condition','mesp_MekMut',"mesp_dnFGFR"),
    c('condition','mesp_MekMut',"control"),
    c('condition','mesp_dnFGFR',"control")
  ),
  test='LRT',reduced=~0+KO,cooksCutoff=F,alpha=.05
)

# 18hpf model, not used
design18 <- subset(b75expDesign,time=='18hpf')
get.dds(
  filt[,row.names(design18)],design18,~0+condition+method,
  path='peakomenocooks/cond18',
  comparisons=list(
    c('condition','handr_MekMut',"handr_dnFGFR"),
    c('condition','handr_MekMut',"control"),
    c('condition','handr_dnFGFR',"control")
  ),
  test='LRT',reduced=~0+method,cooksCutoff=F,alpha=.05
)

# 15-20hpf model
designHandr <- subset(b75expDesign,time%in%c('15hpf','18hpf','20hpf'))
get.dds(
  filt[,row.names(designHandr)],designHandr,~0+condition+time+method,
  path='peakomenocooks/handr',
  comparisons=list(
    c('condition','handr_MekMut',"handr_dnFGFR"),
    c('condition','handr_MekMut',"control"),
    c('condition','handr_dnFGFR',"control")
  ),
  test='LRT',reduced=~0+time+method,cooksCutoff=F,alpha=.05
)

# 10hpf model for only FoxF KO, not used
designFoxf <- subset(b75expDesign,time='10hpf'&condition!='mesp_dnFGFR')
get.dds(
  filt[,row.names(designFoxf)],designFoxf,~0+condition+KO+batch,
  path='peakomenocooks/foxf',
  comparisons=list(
    c('condition','FoxF_KO',"control")
  ),
  test='LRT',reduced=~0+KO+batch,cooksCutoff=F,alpha=.05
)

# time-dependent model without 15 and 20hpf, not used
designLacz <- subset(b75expDesign,condition=='control'&!KO&!time%in%c('15hpf','20hpf'))
get.dds(
  filt[,row.names(designLacz)],designLacz,~0+time+method,
  path='peakomenocooks/lacz6.10.18',
  # pairwise comparisons for LRT
  comparisons=list(
    c('time','6hpf',"10hpf"),
    c('time','10hpf',"18hpf")
  ),
  test='LRT',reduced=~0+method,cooksCutoff=F,alpha=.05
)

# time-dependent model
designLacz <- subset(b75expDesign,condition=='control'&!KO)
get.dds(
  filt[,row.names(designLacz)],designLacz,~0+time+method,
  path='peakomenocooks/lacz',
  # pairwise comparisons for LRT
  comparisons=list(
    c('time','6hpf',"10hpf"),
    c('time','10hpf',"15hpf"),
    c('time','10hpf',"18hpf")
  ),
  test='LRT',reduced=~0+method,cooksCutoff=F,alpha=.05
)
