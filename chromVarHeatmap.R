library(chromVAR)
library(BSgenome.Cintestinalis.KH.KH2013)
library(TFBSTools)
library(BiocParallel)
library(motifmatchr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggseqlogo)

source('data/dirfns.R')
source("data/chromVarPlotFns.R")

writeChromVarHmap("2018-12-12/chromVAR/homer/mesp_all",.01,2)
writeChromVarHmap("2018-12-12/chromVAR/homer/mesp_all",.01,3)

load("2019-01-22/chromVAR/homer/denovoASM/deviations.Rdata")
asmDev <- dev
load("2019-01-22/chromVAR/homer/denovoCardiac/deviations.Rdata")
cardiacDev <- dev

asmDiff <- differentialDeviations(asmDev,'condtime')
cardiacDiff <- differentialDeviations(cardiacDev,'condtime')

asm.z <- deviationScores(asmDev)
colnames(asm.z) <- colData(asmDev)$Name
asm.z <- asm.z[,c(-6,-12:-15)]

cardiac.z <- deviationScores(cardiacDev)
colnames(cardiac.z) <- colData(cardiacDev)$Name
cardiac.z <- cardiac.z[,c(-6,-12:-15)]

asmAvg <- avgZ(asmDev)
cardiacAvg <- avgZ(cardiacDev)

sel <- (asmDiff$p_value_adjusted<.01|
  cardiacDiff$p_value_adjusted<.01)&
  (apply(abs(asmAvg)>1.5,1,any)|
     apply(abs(cardiacAvg)>1.5,1,any))&
  !homer.family%in%c(
    "ERF","AP2EREBP","BBRBPC","C2C2dof","Stat","ZFHD","Myb","MYB",
    "MYBrelated","ND","Trihelix","WRKY","POU,Homeobox",
    "promoter","NA"
  )&
  !is.na(homer.family)

asmHmap <- Heatmap(
  asmAvg[sel,],
  cluster_columns = F,
  split=homer.family[sel],
  col=colorRamp2(c(-3,0,3),c('blue','white','red')),
  row_names_gp = gpar(cex=.5),
  row_title_rot = 0,
  row_title_gp = gpar(cex=.8),
  column_names_gp = gpar(cex=.8),
  show_row_names = F,
  column_title ='De novo ASM'
)
cardiacHmap <- Heatmap(
  cardiacAvg[sel,],
  cluster_columns = F,
  split=homer.family[sel],
  col=colorRamp2(c(-3,0,3),c('blue','white','red')),
  row_names_gp = gpar(cex=.5),
  row_title_rot = 0,
  row_title_gp = gpar(cex=.8),
  column_names_gp = gpar(cex=.8),
  column_title = "De novo cardiac"
)

dir.eps(
  "deNovoChromVarResizeAvgORfdr0.01z2",
  height=sum(sel)*.12+1,
  width=ncol(cardiacAvg)*.5+3
)
draw(asmHmap+cardiacHmap)
dev.off()

denovoCardiacHmap <- chromVarHmap(dev,0.01,2)
denovoAsmHmap <- chromVarHmap(dev,0.01,2)

sapply(
  homer.motifs[sapply(tags(homer.motifs),'[[',"DBID.1")%in%c(
    'MyoD',"MyoG","c-Jun-CRE","Foxf1","Nkx6.1","Smad4","RARg","Pax7","Tbx5","Gata4","Gata6","Otx2","CRX"
  )],
  function(x) dir.gg(
    ggseqlogo(Matrix(x))+ggtitle(tags(x)$DBID.1),
    tags(x)$DBID.1,ext='eps'
  )
)

