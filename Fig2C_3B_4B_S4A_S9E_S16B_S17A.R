#Figs. 2C, 3B, 4B, S4A, S9E, S16B, S17A

source('data/dbHeatmap.R')
source('data/sqlfns.R')
source('data/dirfns.R')
source('data/scatterplotfns.R')
source('data/DESeqFns.R')

library(DBI)

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

bulkGS <- getBulkRNA(con)
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
peakGeneAnnotation <- getAnnotation(con)
gene.names <- peakGeneAnnotation$gene.names

rna <- getRnaDat(con)[c(
  "MA_dnFGFR_LacZ_10hpf","FoxF10hpf_LacZ10hpf",
  "condtime_handrdnfgfr18hpf_handrlacz18hpf",
  "condtime_foxfcamras18hpf_handrlacz18hpf"
)]
atac <- getAtacLib(con,c(
  "condition_mesp_dnFGFR_vs_control","condition_FoxF_KO_vs_control",
  "condition_handr_dnFGFR_vs_control","condition_handr_MekMut_vs_control",
  "tissue_B7.5_vs_mesenchyme","condition_mesp_MekMut_vs_control"
))

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

tmp <- setNames(reshape2::melt(
  mapply(
    function(genes,peaks) setdiff(genes,mergeGenePeak(con,genes,peaks)$GeneID),
    genes=bulkGS[c("MAPK18activated","MAPK18inhibited")],
    peaks=peaksets[c('asmAcc','heartAcc')]
  )
),c("GeneID","gene_set"))
tmp$UniqueNAME <- gene.names[as.character(tmp$GeneID),"UniqueNAME"]
dir.tab(tmp,'MAPK18noDA',row.names=F)

tmp <- rbind(
  cbind(mergeGenePeak2(con,bulkGS$MAPK10activated,peaksets$tvcAcc),gene_set='MAPK10activated'),
  cbind(mergeGenePeak2(con,bulkGS$MAPK10inhibited,peaksets$atmAcc),gene_set='MAPK10inhibited')
) 
tmp <- merge(tmp,rna$MA_dnFGFR_LacZ_10hpf,by.x='GeneID',by.y='row.names')
tmp <- merge(tmp,atac$condition_mesp_dnFGFR_vs_control,by.x='PeakID',by.y='row.names')
tmp <- merge(tmp,atac$condition_mesp_MekMut_vs_control,by.x='PeakID',by.y='row.names')
dir.tab(tmp,"mapk10",row.names=F)

setNames(reshape2::melt(
  mapply(
    function(genes,peaks) setdiff(genes,mergeGenePeak2(con,genes,peaks)$GeneID),
    genes=bulkGS[c("MAPK10activated","MAPK10inhibited")],
    peaks=peaksets[c('tvcAcc','atmAcc')]
  )
),c("GeneID","gene_set"))

# Table S2
tmp <- mergeGenePeak(con,bulkGS$MAPK10activated,peaksets$tvcAcc)
tmp$UniqueNAME <- gene.names[tmp$GeneID,"UniqueNAME"]
tmp <- merge(tmp,dbReadTable(con,'peakfeature')[,1:4])
dir.tab(tmp,'mapk10activated_tvcAcc',row.names=F)

quant.lfc <- c(
  quantile(atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,.0075),
  quantile(atac$condition_mesp_dnFGFR_vs_control$log2FoldChange,.9925)
)

mapk10quant <- lapply(
  bulkGS[c(
    "MAPK10activated","MAPK10inhibited"
  )],
  function(x) unique(mergeGenePeak(
    con,
    x,
    # intersect(x,bulkGS$downreg6hpf),
    row.names(sig.sub(atac$condition_mesp_dnFGFR_vs_control,quant.lfc))
  )$GeneID)
)

tmp <- apply(motifMatches(matches)[,intersect(cisbpDat[cisbpDat$GeneName=="FOXF1/2",1],colnames(motifMatches(matches)))],1,any)
tmp <- names(tmp)[tmp]
tmp <- row.names(homer.matches)[motifMatches(homer.matches)[,"Foxf1"]]

# Fig. 2C
dbHeatmap(
  con,rna[1],
  atac[c(
    "condition_mesp_dnFGFR_vs_control",
    "condition_mesp_MekMut_vs_control"
  )],
  mapk10quant[1],
  append(peaksets[c("atmAcc","tvcAcc")],list(
    intersect(tmp,Reduce(union,peaksets[c("atmAcc","tvcAcc")]))
    # tmp
  )),
  c('gray28','forestgreen','firebrick'),
  list(scrna=prime.denovo),list(scrna=cols),
  'mapk10tvc',
  peak.pch = c(18,18,1),
  peak.cex = c(1.2,1.2,1.2)
)

# Fig. S4A
dbHeatmap(
  con,rna[1],
  atac[c(
    "condition_mesp_dnFGFR_vs_control",
    "condition_mesp_MekMut_vs_control"
  )],
  mapk10quant[2],
  append(peaksets[c("atmAcc","tvcAcc")],list(
    intersect(tmp,Reduce(union,peaksets[c("atmAcc","tvcAcc")]))
    # tmp
  )),
  c('gray28','forestgreen','firebrick'),
  list(scrna=prime.denovo),list(scrna=cols),
  'mapk10atm',
  peak.pch=c(18,18,1)
)

# Fig. 3D, S9E
dbHeatmap(
  con,rna[1:2],atac["condition_FoxF_KO_vs_control"],
  bulkGS[c(
    "FoxFactivated","FoxFinhibited"
  )],
  list(
    intersect(peaksets$closedFoxf,peaksets$tvcAcc),
    # intersect(intersect(peaksets$closedFoxf,peaksets$tvcAcc),tmp)
    intersect(
      intersect(peaksets$closedFoxf,peaksets$tvcAcc),
      row.names(selex.matches)[motifMatches(selex.matches)[,"FoxF"]]
    )
  ),
  c('forestgreen','firebrick'),
  list(scrna=prime.denovo),list(scrna=cols),
  'foxf',
  gene.peak.intersect = F,
  peak.pch = c(18,1)
)

# Fig. 4B
dbHeatmap(
  con,rna,atac[3:4],
  lapply(
    prime.denovo[c('denovoCardiac')],
    intersect,
    union(
      row.names(sig.sub(rna$condtime_foxfcamras18hpf_handrlacz18hpf)),
      row.names(sig.sub(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf))
    )
  ),
  peaksets[c("heartAcc","asmAcc")],c('red','blue'),
  list(scrna=prime.denovo),list(scrna=cols),
  'mapk18denovoCardiac',
  gene.peak.intersect = F
)
# Fig. S17A 
dbHeatmap(
  con,rna,atac[3:4],
  lapply(
    prime.denovo['denovoASM'],
    intersect,
    union(
      row.names(sig.sub(rna$condtime_foxfcamras18hpf_handrlacz18hpf)),
      row.names(sig.sub(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf))
    )
  ),
  peaksets[c("heartAcc","asmAcc")],c('red','blue'),
  list(scrna=prime.denovo),list(scrna=cols),
  'mapk18denovoAsm',
  gene.peak.intersect = F
)

# Fig. S16B
dbHeatmap(
  con,rna,atac[3:4],
  list(
    MAPKinhibited=union(bulkGS$MAPK18inhibited[
      order(rna$condtime_foxfcamras18hpf_handrlacz18hpf[
        bulkGS$MAPK18inhibited,'log2FoldChange'
      ])[1:50]
    ],intersect(prime.denovo$TVCP,bulkGS$MAPK18inhibited)),
    MAPKactivated=union(bulkGS$MAPK18activated[
      order(rna$condtime_handrdnfgfr18hpf_handrlacz18hpf[
        bulkGS$MAPK18activated,'log2FoldChange'
      ])[1:50]
    ],intersect(prime.denovo$STVC,bulkGS$MAPK18activated))
  ),
  peaksets[c("heartAcc","asmAcc")],
  c('red','blue'),
  list(scrna=prime.denovo),list(scrna=cols),
  'mapk18top50',
  gene.peak.intersect = F
)
