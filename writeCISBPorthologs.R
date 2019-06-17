library(DBI)
library(biomaRt)

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

ensembl <- read.delim('KH-ENS.blast',stringsAsFactors = F,header=F)
cisbp.geneid <- read.delim('CISBP/Ciona_intestinalis_2017_11_01_3_41_pm/TF_Information_all_motifs_plus.txt',stringsAsFactors = F)
mart <- useMart('ensembl','cintestinalis_gene_ensembl')
# listAttributes(mart)
transcriptToGene <- select(mart,ensembl[,1],c('ensembl_gene_id','ensembl_transcript_id'),"ensembl_transcript_id")
ensembl.gene <- data.frame(ensembl_transcript_id=ensembl[,1],GeneID=sub('.v.*','',ensembl[,2]))
ensembl.gene <- merge(ensembl.gene,transcriptToGene,'ensembl_transcript_id')
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene[,-1]),-1]
cisbp.geneid <- merge(cisbp.geneid,ensembl.gene,by.x='DBID',by.y='ensembl_gene_id')
cisbp.geneid$GeneID <- sub("KH2012","KH2013",cisbp.geneid$GeneID)
cisbp.geneid$GeneName <- gene.names[cisbp.geneid$GeneID,]
dir.tab(cisbp.geneid,'cisbpGeneID')

cisbpDat <- cisbp.geneid[,c("Motif_ID","GeneID","GeneName","Family_Name","Motif_Type")]
cisbpDat <- cisbpDat[!duplicated(cisbpDat)&cisbpDat$Motif_ID!='.',]
cisbpDat <- rbind(cisbpDat,c("EBF","KH2013:KH.L24.10", "EBF1/2/3/4","EBF","HOMER"))

dbWriteTable(con,"cisbp_orthologs",cisbpDat)
