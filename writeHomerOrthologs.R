library(TFBSTools)
library(DBI)

source("data/dirfns.R")
source('data/sqlfns.R')
source('data/chromVarFns.R')
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ann <- getAnnotation(con)
known.motifs <- getHomerMotifs('known.motifs')
tf.family <- sapply(tags(known.motifs),'[[',"Family_Name")

gene.names <- ann$gene.names
gene.names$regex <- sub('.*_','',gene.names$UniqueNAME)
gene.names$regex <- gsub(
  '([A-Z])|','\\1',
  gsub(
    '/','|',
    gsub(
      "([0-9/]+)","\\(\\1\\)",gene.names$regex
    )))
gene.names['KH2013:KH.C14.307','regex'] <- "MYF5|MYOD"
gene.names['KH2013:KH.L5.5','regex'] <- "ATF3|JDP2"
gene.names['KH2013:KH.C3.17','regex'] <- "RORB|RORGT"
gene.names['KH2013:KH.L24.10','regex'] <- "EBF"
gene.names[grep("CREB",gene.names$UniqueNAME),"regex"] <- "CRE"
gene.names["KH2013:KH.C7.127","regex"] <- "FOXP"
gene.names["KH2013:KH.C1.1116",'regex'] <- "HAND"
gene.names[grep("FOXL",gene.names$UniqueNAME),"regex"] <- "FOXL"
# gene.names["KH2013:KH.C14.604",'family'] <- "bHLH"
# gene.names["KH2013:KH.L152.2",'family'] <- "Homeobox"
# gene.names["KH2013:KH.C1.1116",c('UniqueNAME','family','regex')] <- c("HAND-R","bHLH","HAND")

# sel <- grep("FOX",gene.names$UniqueNAME)
# gene.names[sel,"regex"] <- sub('\\(.*','[0-9]+',gene.names[sel,"regex"])
gene.names$regex <- paste0("(^|[[:space:]]|[[:punct:]]|-)",gene.names$regex,"($|[[:space:]]|[[:punct:]]|-)")

gene.motif <- sapply(gene.names$regex,grep,names(known.motifs),ignore.case=T)
sel <- sapply(gene.motif,length)>0
gene.motif <- do.call(rbind,mapply(
  cbind,
  lapply(gene.motif[sel],function(x) cbind(
    ID=names(known.motifs)[x],
    Family_Name=tf.family[x]
  )),
  GeneID=row.names(gene.names)[sel],
  GeneName=gene.names$UniqueNAME[sel]
))

ma <- dbReadTable(con,'microarray')
ma <- subset(ma,probeannotation11__kh!='')
row.names(ma) <- paste0('KH2013:',ma$probeannotation11__kh)

motif.gene <- sapply(
  paste0('(_|,|\\s|^)',gsub('\\.','-',names(known.motifs)),'(_|,|\\s|$)'),
  grep,
  ma$NewAnnotation_ASMTVC2_NAME_Best_Hit_Human_Proteome,
  ignore.case=T
)
motif.gene <- lapply(motif.gene,function(x) ma[
  x,c("GeneID","GeneName","NewAnnotation_ASMTVC2_NAME_Best_Hit_Human_Proteome")
  ])
sel <- sapply(motif.gene,nrow)>0
motif.gene <- mapply(
  cbind,
  ID=names(known.motifs)[sel],
  Family_Name=tf.family[sel],
  motif.gene[sel],
  SIMPLIFY = F
)
motif.gene <- do.call(rbind,motif.gene)
dir.tab(motif.gene[motif.gene$ID!="Unknown",],'homerToKH')

homerDat <- rbind(gene.motif,motif.gene[,-5])
homerDat <- homerDat[
  !duplicated(homerDat[,c("ID","GeneID")])&homerDat$ID!="Unknown"&!grepl("SeqBias",homerDat$ID),
]

dir.tab(homerDat,'khToHomer')
dbWriteTable(con,'homer_orthologs',homerDat,overwrite=T)
