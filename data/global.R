library(DBI)
library(shiny)
library(DT)

source('sqlfns.R')
source('dirfns.R')
source('dbDiamond.R')
source('DESeqFns.R')

con <- dbConnect(RSQLite::SQLite(),'atacCiona.db')
rna.atac.res <- dbReadTable(con,'rnaseq_atac_lib')

peaksets <- getPeaksets(con)
scrna <- getScRNA(con)
bulkGS <- getBulkRNA(con)
atac <- getAtac(con)
rna <- getRnaDat(con)

prime.denovo <- scrna[c(
  'primedCardiac','primedASM','denovoCardiac','denovoASM','TVCP','STVC','ATM','mesenchyme'
)]
ebf <- scrna[c('ebfActivated','ebfInhibited')]
cols <- c(
  primedCardiac='red',primedASM='blue',
  denovoCardiac='hotpink',denovoASM='skyblue',
  TVCP='forestgreen',STVC='orange',ATM='gray28',
  ebfActivated='lightgreen',ebfInhibited='purple'
)
