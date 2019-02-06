# Figs. 2A, S20

library(chromVAR)
library(BSgenome.Cintestinalis.KH.KH2013)
library(TFBSTools)
library(BiocParallel)
library(motifmatchr)
library(ComplexHeatmap)
library(circlize)

source('data/dirfns.R')
source("data/chromVarPlotFns.R")

con <- dbConnect(RSQLite::SQLite(),"data/atacCiona.db")

# Fig. 2A
mapk10 <- dbReadTable(con,'mapk10ChromVAR')[,-8]
writeDevHmap(mapk10,"")

writeChromVarHmap("2019-01-22/chromVAR/homer/mesp_all",.01,3)
writeChromVarHmap("../hpcscripts/2019-02-01/chromVAR/homer/mesp")

# Fig. S20

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
cardiac.z <- cardiac.z[,c(-6,-8,-9,-12:-15)]

asmAvg <- avgZ(asmDev[,c(-6,-8,-9,-12:-15)])
cardiacAvg <- avgZ(cardiacDev[,c(-6,-8,-9,-12:-15)])

sel <- (asmDiff$p_value_adjusted<.01|
  cardiacDiff$p_value_adjusted<.01)&
  (apply(abs(asmAvg)>1.5,1,any)|
     apply(abs(cardiacAvg)>1.5,1,any))&
  mcols(asmDev)$Family_Name%in%c(
    "bHLH","bZIP","EBF","ETS","Forkhead","Homeobox","MAD","MADS","NR","Paired,Homeobox","T-box","Zf"
  )&
  !is.na(mcols(asmDev)$Family_Name)

asmHmap <- Heatmap(
  asmAvg[sel,],
  cluster_columns = F,
  split=unlist(mcols(asmDev)$Family_Name)[sel],
  col=colorRamp2(c(-3,0,3),c('blue','white','red')),
  row_names_gp = gpar(cex=.5),
  row_title_rot = 90,
  row_title_gp = gpar(cex=.8),
  column_names_gp = gpar(cex=.8),
  show_row_names = F,
  column_title ='De novo ASM'
)
cardiacHmap <- Heatmap(
  cardiacAvg[sel,],
  cluster_columns = F,
  split=unlist(mcols(cardiacDev)$Family_Name)[sel],
  col=colorRamp2(c(-3,0,3),c('blue','white','red')),
  row_names_gp = gpar(cex=.5),
  row_title_rot = 90,
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
