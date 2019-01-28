dbWriteRownamesAs <- function(conn,name,value,row.names='row_names',...){
  require(DBI)
  value <- as.data.frame(value)
  colnames <- c(row.names,names(value))
  value <- do.call(data.frame,append(list(row.names(value)),value))
  names(value) <- colnames
  dbWriteTable(conn,name,value,...,row.names=F)
}

getTable <- function(sql) sub(".*?FROM (.*?) .*","\\1",sql)

getSig <- function(con,table,lfc,p=0.05) {
  require(DBI)
  dbGetQuery(con,paste0(
    "SELECT * FROM ",table," WHERE log2FoldChange",lfc," AND padj<",as.character(p)
  ))
}

getAbsSig <- function(con,table,lfc,p=0.05) {
  require(DBI)
  dbGetQuery(con,paste0(
    "SELECT * FROM ",table," WHERE abs(log2FoldChange)>",as.character(lfc)," AND padj<",as.character(p)
  ))
}

geneToPeak <- function(con,genes){
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT GeneID, PeakID FROM geneToPeak WHERE GeneID IN ("',
    paste(genes,collapse = '", "'),
    '")'
  )))
  # sel <- dbSendQuery(con, 'SELECT PeakID FROM geneToPeak WHERE GeneID=:1',data=genes)
  # return(sel)
  # dbBind(sel,param = list(genes=genes))
  # peaks <- dbFetch(sel)
  # dbClearResult(sel)
  # return(unique(peaks$PeakID))
}

splitByDistinct <- function(con,table,field,select='*',...) {
  require(DBI)
  sapply(
    dbGetQuery(con,paste(
      'SELECT DISTINCT', field, "FROM", table
    ))[,1],
    function(x) dbGetQuery(con,paste0(
      'SELECT ',select,' FROM ',table,' WHERE ',field,'="',x,'"'
    ),...),
    simplify = F
  )
}

peakToGene <- function(con,peaks){
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT PeakID, GeneID FROM geneToPeak WHERE PeakID IN("',
    paste(peaks,collapse = '","'),
    '")'
  )))
  # sel <- dbGetQuery(con, 'SELECT GeneID, PeakID FROM geneToPeak WHERE PeakID=:1',bind.data=peaks)
  # return(sel)
  # sel <- dbSendQuery(con, 'SELECT GeneID FROM geneToPeak WHERE PeakID=:peaks')
  # dbBind(sel,param = list(peaks=peaks))
  # genes <- dbFetch(sel)
  # dbClearResult(sel)
  # return(unique(genes$GeneID))
}

mergeGenePeak <- function(con,genes,peaks){
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT GeneID,PeakID FROM geneToPeak WHERE GeneID IN("',
    paste(genes,collapse = '","'),
    '") AND PeakID IN("',
    paste(peaks,collapse='","'),'")'
  )))
}

mergeGenePeak2 <- function(con,genes,peaks){
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT geneToPeak.GeneID,UniqueNAME,geneToPeak.PeakID,peakfeature.chr,peakfeature.start,peakfeature.end FROM (geneToPeak LEFT JOIN gene_name ON geneToPeak.GeneID=gene_name.GeneID) LEFT JOIN peakfeature ON geneToPeak.PeakID=peakfeature.PeakID WHERE geneToPeak.GeneID IN("',
    paste(genes,collapse = '","'),
    '") AND geneToPeak.PeakID IN("',
    paste(peaks,collapse='","'),
    '")'
  )))
}


# getGeneNames <- function(con) dbReadTable(con,'gene_name',row.names="GeneID")

# may improve speed by reading all association into R environment
getAnnotation <- function(con){
  require(DBI)
  require(GenomicRanges)
  require(BSgenome.Cintestinalis.KH.KH2013)
  dat <- dbReadTable(con,'geneToPeak')
  geneToPeak <- split(dat$PeakID,dat$GeneID)
  peakToGene <- split(dat$GeneID,dat$PeakID)
  features <- dbReadTable(con,"peakfeature",row.names="PeakID")
  peaks <- GRanges(features[,1],IRanges(features[,2],features[,3],names = row.names(features)))
  features <- features[,-1:-3]
  sapply(1:length(features),function(x) features[,x] <<- as.logical(features[,x]))
  genes <- dbReadTable(con,"gene_name")
  gene.names <- data.frame(
    UniqueNAME=genes$UniqueNAME,row.names = genes$GeneID,stringsAsFactors = F
  )
  genes <- GRanges(
    genes$chr,
    IRanges(genes$start,genes$end,names = genes$GeneID),
    UniqueNAME=genes$UniqueNAME
  )
  genewindow <- resize(genes,width(genes)+20000,fix = "center")
  # trim out-of-bounds ranges
  start(genewindow)[start(genewindow) < 1] <- 1
  # trim ranges past the chromosome length
  mapply(
    function(chr,chrlen){
      end(genewindow)[
        which(seqnames(genewindow)==chr&end(genewindow)>chrlen)
      ] <<- chrlen
    },
    seqnames(BSgenome.Cintestinalis.KH.KH2013),
    seqlengths(BSgenome.Cintestinalis.KH.KH2013)
  )
  return(list(
    geneToPeak=geneToPeak,
    peakToGene=peakToGene,
    features=features,
    peaks=peaks,
    genes=genes,
    genewindow=genewindow,
    gene.names=gene.names
  ))
}

getScRNA <- function(con){
  require(DBI)
  scrna <- Reduce(append,sapply(c(
    'scrna','ebfdat'
  ),function(x) {
    res <- dbReadTable(con,x)
    return(split(res$GeneID,res$geneset))
  }))
  scrna$ATM <- dbReadTable(con,"ATM_genes_from_ANISEED")$GeneID
  scrna$mesenchyme <- dbReadTable(con,"mesenchyme20hpf")$GeneID
  scrna$Cardiac <- dbGetQuery(
    con,'SELECT DISTINCT GeneID FROM scrna WHERE geneset IN ("Pancardiac","FHP","FHP14","SHP")'
  )[,1]
  scrna$primedCardiac <- dbGetQuery(
    con,'SELECT DISTINCT GeneID FROM scrna WHERE geneset IN ("Pancardiac","FHP","FHP14","SHP") AND Type="Primed"'
  )[,1]
  scrna$primedASM <- dbGetQuery(
    con,'SELECT DISTINCT GeneID FROM scrna WHERE geneset="ASM" AND Type="Primed"'
  )[,1]
  scrna$denovoCardiac <- dbGetQuery(
    con,'SELECT DISTINCT GeneID FROM scrna WHERE geneset IN ("Pancardiac","FHP","FHP14","SHP") AND Type="De Novo"'
  )[,1]
  scrna$denovoASM <- dbGetQuery(
    con,'SELECT DISTINCT GeneID FROM scrna WHERE geneset="ASM" AND Type="De Novo"'
  )[,1]
  scrna$STVC <- setdiff(scrna$STVC,scrna$TVCP)
  scrna$ATM <- setdiff(scrna$ATM,Reduce(union,scrna[c(
    "FHP","FHP14","SHP","Pancardiac","STVC","TVCP"
  )]))
  scrna$mesenchyme <- setdiff(scrna$mesenchyme,Reduce(union,scrna[c(
    "ASM","FHP","FHP14","SHP","Pancardiac","STVC","TVCP","ATM"
  )]))
  return(scrna)
}

getBulkRNA <- function(con){
  require(DBI)
  bulkGS <- dbReadTable(con,'bulkRNAgenesets',row.names='geneset')
  bulkGS <- apply(bulkGS,1,function(x) dbGetQuery(con,x)$GeneID)
  return(bulkGS)
}

getRnaDat <- function(con){
  require(DBI)
  res <- splitByDistinct(
    con,'handr_rnaseq','comparison','GeneID,log2FoldChange,padj',row.names='GeneID'
  )
  # foxf <- splitByDistinct(
  #   con,'foxf_rnaseq','comparison','GeneID,log2FoldChange,padj',row.names='GeneID'
  # )
  res$MA_dnFGFR_LacZ_10hpf <- dbGetQuery(
    con,
    'SELECT x AS GeneID,microarray.Mesp_dnfgfrvwt_logfc AS log2FoldChange,microarray.dnfgfrvwt_pval AS padj FROM (SELECT DISTINCT GeneID AS x FROM handr_rnaseq) LEFT JOIN microarray ON x=microarray.GeneID',
    row.names="GeneID"
  )
  res$FoxF10hpf_LacZ10hpf <- dbGetQuery(
    con,
    'SELECT x AS GeneID,log2FoldChange,padj FROM (SELECT DISTINCT GeneID AS x FROM handr_rnaseq) LEFT JOIN (SELECT * FROM foxf_rnaseq WHERE comparison="FoxF10hpf_LacZ10hpf") ON x=GeneID',
    row.names="GeneID"
  )
  res$MA_6hpf_10hpf <- dbGetQuery(
    con,
    'SELECT x AS GeneID,microarray.lgvwt_logfc AS log2FoldChange,microarray.lgvwt_pval AS padj FROM (SELECT DISTINCT GeneID AS x FROM handr_rnaseq) LEFT JOIN microarray ON x=microarray.GeneID',
    row.names="GeneID"
  )
  return(res)
}

getPeaksets <- function(con){
  require(DBI)
  peaksets <- dbReadTable(con,'peaksets',row.names='peakset')
  peaksets <- apply(peaksets,1,function(x) dbGetQuery(con,x)$PeakID)
  return(peaksets)
}

getAtac <- function(con)  splitByDistinct(
  con,'atacseq','comparison','PeakID,log2FoldChange,padj',row.names='PeakID'
)

getAtacLib <- function(con,lib) {
  require(DBI)
  sapply(lib,function(x) dbGetQuery(
    con,
    paste0('SELECT PeakID,log2FoldChange,padj FROM atacseq WHERE comparison="',x,'"'),
    row.names="PeakID"
  ),simplify = F)
  # dbBind(query,list(lib=lib))
  # return(dbFetch(query,row.names="PeakID"))
}
