melt.rename <- function(dat,value='value',...){
# appends list of data.frames as in reshape2::melt
  res <- sapply(dat,as.data.frame,simplify = F)
  res <- mapply(
    cbind,dat,names(dat),stringsAsFactors=F,SIMPLIFY = F
  )
  res <- lapply(res, function(x) setNames(x,c(names(x)[-length(x)],value)))
  res <- do.call(rbind,res)
  return(res)
}


dbWriteKey <- function(conn,name,value,primary=F,foreign=F,row.names=F,overwrite=T,...){
  # wrapper function for dbWriteTable which allows specification of primary and foreign keys
  require(RSQLite)
  value <- as.data.frame(value)
  name <- make.db.names(conn,name)
  names(value) <- make.db.names(conn,names(value))
  if(is.character(row.names)){
    colnames <- c(row.names,names(value))
    value <- do.call(data.frame,append(list(row.names(value)),value))
    names(value) <- colnames
  }
  s <- sprintf("CREATE TABLE %s(%s", name, paste(names(value), collapse = ", "))
  if(is.character(primary)) s <- paste0(s,', PRIMARY KEY(',primary,')')
  if(is.character(foreign)) s <- paste0(s,paste0(
    ', FOREIGN KEY(',
    sub('.*\\((.*)\\)',"\\1",foreign),
    ') REFERENCES ',
    foreign,
    collapse = ''
    ))
  s <- paste0(s,')')
  if(overwrite&dbExistsTable(conn,name)) dbSendStatement(conn,paste("DROP TABLE",name))
  dbSendStatement(conn, s)
  dbWriteTable(conn,name,value,...,row.names=F,append=T)
}

dbWriteGenes <- function(...) dbWriteKey(...,foreign = 'gene_name(GeneID)')
dbWritePeaks <- function(...) dbWriteKey(...,foreign = 'peakfeature(PeakID)')

dbWriteRownamesAs <- function(conn,name,value,row.names='row_names',...){
  # wrapper function for dbWriteTable allowing a user-specified name for row.names attriute
  require(DBI)
  value <- as.data.frame(value)
  colnames <- c(row.names,names(value))
  value <- do.call(data.frame,append(list(row.names(value)),value))
  names(value) <- colnames
  dbWriteTable(conn,name,value,...,row.names=F)
}

# extracts table name from SQL query
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
  # accepts a vector of GeneIDs
  # returns a data.frame of all peaks associated to genes
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
  # accepts a database connection, table name and field name
  # splits table by field into a list of data.frames
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
  # accepts a database connection and a vector of PeakIDs
  # returns a data.frame of all GeneIDs associated to peaks
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
  # accepts a database connection, list of GeneIDs, and list of PeakIDs
  # returns a data.frame of all associations between peaks and genes
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT GeneID,PeakID FROM geneToPeak WHERE GeneID IN("',
    paste(genes,collapse = '","'),
    '") AND PeakID IN("',
    paste(peaks,collapse='","'),'")'
  )))
}

mergeGenePeak2 <- function(con,genes,peaks){
  # as mergeGenePeak, but also returns gene names and peak coordinates
  require(DBI)
  return(dbGetQuery(con,paste0(
    'SELECT geneToPeak.GeneID,UniqueNAME,geneToPeak.PeakID,peakfeature.chr,peakfeature.start,peakfeature.end FROM (geneToPeak LEFT JOIN gene_name ON geneToPeak.GeneID=gene_name.GeneID) LEFT JOIN peakfeature ON geneToPeak.PeakID=peakfeature.PeakID WHERE geneToPeak.GeneID IN("',
    paste(genes,collapse = '","'),
    '") AND geneToPeak.PeakID IN("',
    paste(peaks,collapse='","'),
    '")'
  )))
}

# returns a 1-column dataframe of gene names with GeneID as row.names
getGeneNames <- function(con) dbGetQuery(con,'SELECT GeneID,UniqueNAME FROM gene_name',row.names="GeneID")

getAnnotation <- function(con){
# it may improve speed to read all association into R environment
  # reads all gene-peak associations from a database connection.
  # returns a list with attributes:
    # geneToPeak  list of vectors of PeakIDs split by GeneID
    # peakToGene  list of vectors of GeneID split by PeakID
    # features  logical matrix of genomic features overlapping peaks
    # peaks   GRanges of peakome
    # genes   GRanges of KH2013 gene loci
    # genewindow GRanges of KH2013 gene loci +/- 10kbp from TSS/TTS
    # gene.names  mapping of GeneID to gene name
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
  # accepts a database connection
  # returns a list of GeneID vectors split by gene sets from scRNAseq, 
  # EBF perturbation microarray, ATM genes from ANISEED, and mesenchymal genes
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
  # accepts a database connection
  # returns a list of GeneID vectors split by DE gene sets from bulk RNA-seq or microarrays
  # the selection of DE genes for each set is given in the table 'bulkRNAgenesets'
  require(DBI)
  bulkGS <- dbReadTable(con,'bulkRNAgenesets',row.names='geneset')
  bulkGS <- apply(bulkGS,1,function(x) dbGetQuery(con,x)$GeneID)
  return(bulkGS)
}

getRnaDat <- function(con){
  # accepts a database connection
  # returns a list of data.frames, each of which is a pairiwise comparison from bulk RNA-seq or microarrays
  require(DBI)
  res <- splitByDistinct(
    con,'handr_rnaseq','comparison','GeneID,log2FoldChange,padj',row.names='GeneID'
  )
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
  # accepts a database connection
  # returns a list of PeakID vectors split by peak sets
  # the selection of peaks in each peak set is given in the 'peaksets' table
  require(DBI)
  peaksets <- dbReadTable(con,'peaksets',row.names='peakset')
  peaksets <- apply(peaksets,1,function(x) dbGetQuery(con,x)$PeakID)
  return(peaksets)
}

getAtac <- function(con)  splitByDistinct(
  # accepts a database connection
  # returns a list of data.frames, each of which is a pairwise comparison from ATAC-seq
  con,'atacseq','comparison','PeakID,log2FoldChange,padj',row.names='PeakID'
)

getAtacLib <- function(con,lib) {
  # accepts a database connection and the name of an ATAC-seq comparison in the table atacseq
  # returns a data.frame with the ATAC-seq comparison specified
  require(DBI)
  sapply(lib,function(x) dbGetQuery(
    con,
    paste0('SELECT PeakID,log2FoldChange,padj FROM atacseq WHERE comparison="',x,'"'),
    row.names="PeakID"
  ),simplify = F)
}
