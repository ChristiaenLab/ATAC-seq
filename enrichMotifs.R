library(TFBSTools)
library(motifmatchr)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)
library(chromVAR)
library(SummarizedExperiment)

source('data/chromVarFns.R')
source('data/sqlfns.R')
source('data/getMotifs.R')
source('data/alignMotifs.R')
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ann <- getAnnotation(con)
peaksets <- getPeaksets(con)

bg <- letterFrequency(Views(Cintestinalis,ann$peaks),c("A","C","G","T"))
bg <- apply(bg,2,sum)
bg <- bg/sum(bg)

motifs <- reduceMotifs(con,F,F)
motif.dup <- setNames(motifs,name(motifs))
alignMotifs('nkx2_3.fa',motif.dup,bg=bg)
alignMotifs('smurf.fa',motif.dup,bg=bg)
alignMotifs('hand.fa',motif.dup,bg=bg)
alignMotifs('handfull.fa',motif.dup,bg=bg)

motifs <- reduceMotifs(con,F,F,F)
motif.dup <- setNames(motifs,name(motifs))
alignMotifs('nkx2_3.fa',suffix = 'dup',motif.dup,bg=bg)
alignMotifs('smurf.fa',suffix = 'dup',motif.dup,bg=bg)
alignMotifs('hand.fa',suffix = 'dup',motif.dup,bg=bg)
alignMotifs('handfull.fa',suffix = 'dup',motif.dup,bg=bg)

matches <- matchMotifs(motifs,ann$peaks,Cintestinalis,bg=bg)
sel <- split(names(motifs),ID(motifs))

fn <- function(x){
  if(length(x)>1){
    pwms <- combn(x,2)
    pearson <- apply(pwms,2,function(x) PWMSimilarity(
      motifs[[x[1]]],
      motifs[[x[2]]],
      "Pearson"
    ))
    return(x[!x%in%pwms[2,pearson>.90]])
  } else return(x)
}

motif.cor <- lapply(sel, fn)

dupct <- lapply(sel,function(x) apply(motifMatches(matches)[,x,drop=F],2,sum))
sel <- sapply(dupct,function(x) names(x)[which.max(x)])

sel <- unlist(motif.cor)
tf.name <- sapply(tags(motifs[sel]),'[[',"DBID.1")
tf.family <- sapply(tags(motifs[sel]),'[[',"Family_Name")
tf.kh <- ID(motifs[sel])

load('2019-06-25/chromVarOut.Rdata')

mespDev <- mespDev[,c(-6,-12:-15)]
diff_acc <- differentialDeviations(mespDev,'condtime')

sel <- split(names(motifs),ID(motifs))
sel <- sapply(sel,function(x) x[which.min(diff_acc[x,1])])

dev <- deviationScores(mespDev)
colnames(dev) <- colData(mespDev)$Name
dev <- dev[
  # diff_acc$p_value_adjusted<.05&apply(dev,1,function(x) any(abs(x)>1.5)),
  sel,
]

# motifs <- reduceMotifs(con)
tf.name <- sapply(tags(motifs[sel]),'[[',"DBID.1")
tf.family <- sapply(tags(motifs[sel]),'[[',"Family_Name")
# matches <- matchMotifs(motifs,ann$peaks,Cintestinalis,bg=bg)

hyper <- lapply(append(peaksets[c(5,6,1,2)],setNames(
  lapply(peaksets[c('open6','closed6')],intersect,peaksets$tvcAcc),
  c("earlyTVC","lateTVC")
),3),lHyper,motifMatches(matches)[,sel])

odds <- sapply(hyper,'[',,'log2OddsRatio')
fdr <- sapply(hyper,'[',,'padj')
odds[fdr>.05] <- 0
# row.names(odds) <- names(motifs)
odds <- cbind(TFGeneID=ID(motifs[sel]),TF=names(motifs[sel]),family=tf.family,as.data.frame(odds),stringsAsFactors=F)
odds <- odds[!apply(odds[,-1:-3]<1.5,1,all),c(rep(T,3),!apply(odds[,-1:-3]<1.5,2,all))]

mat <- merge(odds,dev,'row.names')[,-1]

tf.kh.gene <- strsplit(mat$TFGeneID,';')
times <- sapply(tf.kh.gene,length)
mat <- mapply(
  function(x,y,z) cbind(
    GeneID=z,
    do.call(rbind,replicate(x,y,F)),
    stringsAsFactors=F
  ),
  times,
  split(mat,factor(mat$TF,mat$TF)),
  tf.kh.gene,
  SIMPLIFY = F
)
# mat <- mapply(cbind,mat,GeneID=tf.kh.gene,stringsAsFactors=F,SIMPLIFY = F)
# mat <- do.call(rbind,unlist(mat,F))
# mat <- Reduce(function(x,y) {
#   x[tf.kh.gene[[y]]] <- mat[[y]]
#   return(x)
# },1:length(mat),list())
mat <- do.call(rbind,mat)


ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
names(ma) <- as.character(seq(6,20,2))
mat <- merge(mat,ma,by.x="GeneID",by.y="row.names")
row.names(mat) <- make.unique(as.character(mat$TF))
# row.names(mat) <- gene.names[mat$GeneID,]

# heatmap_width <- .25+max(nchar(row.names(mat)))*.05+max(nchar(as.character(mat$family)))*.05+(ncol(mat)-3)/4+heatmap_width
heatmap_height <- nrow(mat)*.20+max(nchar(names(mat)))*.08
lwidth <- max(nchar(mat$family))*.08
rwidth <- max(nchar(row.names(mat)))*.08

hm2 <- Heatmap(
  as.matrix(mat[,names(ma)]),
  split=mat$family,
  heatmap_width = unit(ncol(ma)*.20+.5,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  cluster_columns = F,
  cluster_rows = T,
  cluster_row_slices = T,
  show_row_dend = F,
  row_dend_side='right',
  row_title_rot = 0,
  row_title_gp = gpar(cex=1),
  column_gap=unit(.2,'in'),
  name="TF expression \n(log2FC)"
)

hm1 <- Heatmap(
  as.matrix(mat[,names(odds)[c(-1:-3)]]),
  # breaks = c(0,4),
  split=mat$family,
  cluster_columns = F,
  cluster_row_slices=T,
  show_column_dend=F,
  show_row_dend=F,
  row_dend_side='right',
  heatmap_width = unit(rwidth+(ncol(odds)-3)*.20+.5,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  # heatmap_width = unit((ncol(odds)-2)/4+heatmap_width,'in'),
  # heatmap_height = unit(heatmap_height,'in')
  row_title_rot = 0,
  row_title_gp = gpar(cex=1),
  col=colorRamp2(c(0,3),c('white','black')),
  column_gap=unit(.2,'in'),
  name='TF motif enrichment \n(log2OR)'
)

hm3 <- Heatmap(
  as.matrix(mat[,colnames(dev)]),
  # breaks = c(0,4),
  split=mat$family,
  cluster_columns = F,
  cluster_row_slices=T,
  show_column_dend=F,
  show_row_dend=F,
  row_dend_side='right',
  heatmap_width = unit(lwidth+(ncol(dev))*.20+.5,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  # heatmap_width = unit((ncol(odds)-2)/4+heatmap_width,'in'),
  # heatmap_height = unit(heatmap_height,'in')
  row_title_rot = 0,
  row_title_gp = gpar(cex=1),
  col = colorRamp2(c(-5,0,5),c("blue","white","red")),
  column_gap=unit(.2,'in'),
  name="TF motif accessibility \n(deviation z-score)"
)

dir.eps(
  "or1.5corfilt",
  width=(ncol(mat)-3)*.25+2+.25+lwidth+rwidth,
  height=heatmap_height+4
)
draw(hm3+hm2+hm1)
dev.off()

hm3 <- Heatmap(
  dev[row.names(mat),],split=tf.family[row.names(dev)]
)

writeChromVarHmap("2019-06-25/",.01,3)

peakHyper(append(peaksets[1:12],list(
  earlyTVC=intersect(peaksets$tvcAcc,peaksets$open6),
  lateTVC=intersect(peaksets$tvcAcc,peaksets$closed6)
)),'peaksets',motifMatches(matches),tf.family)

row.names(mat) <- gene.names[mat$GeneID,]
hm1 <- Heatmap(
  as.matrix(mat[,names(ma)]),
  split=mat$family,
  heatmap_width = unit((ncol(ma))*.25+
    max(nchar(row.names(mat)))*.05+
    max(nchar(as.character(mat$family)))*.05+4,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  cluster_columns = F,
  cluster_row_slices = F,
  row_title_rot = 0,
  row_title_gp = gpar(cex=1),
  right_annotation = rowAnnotation(TF=row_anno_text(mat$TF))
)

dir.eps(
  "tmp",
  width=(ncol(ma)-3)*.25+6+.25+
    max(nchar(row.names(mat)))*.05+
    max(nchar(as.character(mat$family)))*.05,
  height=heatmap_height+4)
draw(hm1)
dev.off()
