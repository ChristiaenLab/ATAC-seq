library(TFBSTools)
library(motifmatchr)
library(DBI)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(ComplexHeatmap)
library(circlize)

source('data/chromVarFns.R')
source('data/sqlfns.R')
source('data/getMotifs.R')
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ann <- getAnnotation(con)
motifs <- reduceMotifs(con)
tf.name <- sapply(tags(motifs),'[[',"DBID.1")
tf.family <- sapply(tags(motifs),'[[',"Family_Name")

bg <- letterFrequency(Views(Cintestinalis,ann$peaks),c("A","C","G","T"))
bg <- apply(bg,2,sum)
bg <- bg/sum(bg)
matches <- matchMotifs(motifs,ann$peaks,Cintestinalis,bg=bg)

hyper <- lapply(append(peaksets[1:12],setNames(
  lapply(peaksets[c('open6','closed6')],intersect,peaksets$tvcAcc),
  c("earlyTVC","lateTVC")
)),lHyper,motifMatches(matches))

odds <- sapply(hyper,'[',,'log2OddsRatio')
fdr <- sapply(hyper,'[',,'padj')
odds[fdr>.05] <- 0
# row.names(odds) <- names(motifs)
odds <- cbind(TF=names(motifs),family=tf.family,as.data.frame(odds),stringsAsFactors=F)

tf.kh.gene <- strsplit(tf.kh,';')
times <- sapply(tf.kh.gene,length)
odds <- mapply(
  function(x,y,z) cbind(
    GeneID=z,
    do.call(rbind,replicate(x,y,F)),
    stringsAsFactors=F
  ),
  times,
  split(odds,factor(odds$TF,odds$TF)),
  tf.kh.gene,
  SIMPLIFY = F
)
# odds <- mapply(cbind,odds,GeneID=tf.kh.gene,stringsAsFactors=F,SIMPLIFY = F)
# odds <- do.call(rbind,unlist(odds,F))
# odds <- Reduce(function(x,y) {
#   x[tf.kh.gene[[y]]] <- odds[[y]]
#   return(x)
# },1:length(odds),list())
odds <- do.call(rbind,odds)

odds <- odds[!apply(odds[,-1:-3]<1.5,1,all),c(rep(T,3),!apply(odds[,-1:-3]<1.5,2,all))]

ma <- dbReadTable(con,'microarray',row.names="GeneID")[,23:30]
names(ma) <- as.character(seq(6,20,2))
mat <- merge(odds,ma,by.x="GeneID",by.y="row.names")
row.names(mat) <- make.unique(as.character(mat$TF))
# row.names(mat) <- gene.names[mat$GeneID,]

# heatmap_width <- .25+max(nchar(row.names(mat)))*.05+max(nchar(as.character(mat$family)))*.05+(ncol(mat)-3)/4+heatmap_width
heatmap_height <- nrow(mat)/4+max(nchar(names(mat)))*.05
lwidth <- max(nchar(mat$family))*.05
rwidth <- max(nchar(row.names(mat)))*.05

hm2 <- Heatmap(
  as.matrix(mat[,names(ma)]),
  split=mat$family,
  heatmap_width = unit(ncol(ma)*.25+rwidth+.5,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  cluster_columns = F,
  cluster_rows = T,
  cluster_row_slices = T,
  show_row_dend = F,
  row_dend_side='right',
  row_title_rot = 0,
  row_title_gp = gpar(cex=1)
)
hm1 <- abs.hmap(
  as.matrix(mat[,names(odds)[c(-1:-3)]]),
  breaks = c(0,4),
  split=mat$family,
  cluster_row_slices=T,
  show_column_dend=F,
  show_row_dend=F,
  row_dend_side='right',
  heatmap_width = unit(lwidth+(ncol(odds)-3)*.25+.5,'in'),
  heatmap_height = unit(heatmap_height,'in'),
  # heatmap_width = unit((ncol(odds)-2)/4+heatmap_width,'in'),
  # heatmap_height = unit(heatmap_height,'in')
  row_title_rot = 0,
  row_title_gp = gpar(cex=1)
)
dir.eps(
  "tmp2",
  width=(ncol(ma)-3)*.25+6+.25+lwidth+rwidth,
  height=heatmap_height+4
)
draw(hm1+hm2)
dev.off()

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
