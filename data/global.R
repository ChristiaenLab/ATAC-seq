library(DBI)
library(shiny)

source('sqlfns.R')
source('dirfns.R')
source('diamondplotFns.R')
source('dbDiamond.R')

con <- dbConnect(RSQLite::SQLite(),'atacCiona.db')
rna.atac.res <- dbReadTable(con,'rnaseq_atac_lib')
