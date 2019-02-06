library(motifmatchr)
library(TFBSTools)
library(DBI)
library(BSgenome.Cintestinalis.KH.KH2013)

source('data/plotMotifs.R')
source('data/sqlfns.R')
source('data/chromVarFns.R')
source('data/motifHyperFns.R')
source('data/dirfns.R')
source('data/corHeatmap.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')
scrna <- getScRNA(con)
peaksets <- getPeaksets(con)
ann <- getAnnotation(con)

known.motifs <- getHomerMotifs('known.motifs')
matches <- matchMotifs(known.motifs,ann$peaks,BSgenome.Cintestinalis.KH.KH2013,'subject','matches')
tf.family <- sapply(tags(known.motifs),'[[',"Family_Name")

tmp <- stack(sapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(
        Reduce(union,peaksets[c("open6",'closed6','closed18')]),
        peaksets$open18
      ),
      denovoPeaks=peaksets$open18
    ),
    function(y) mergeGenePeak2(con,x,y)
  )
))
dir.tab(tmp[,-3],'primed_denovo_cardiac_asm',col.names=F,row.names=F)

mekmut.dnfgfr.18 <- getAtacLib(con,c('condition_handr_MekMut_vs_control','condition_handr_dnFGFR_vs_control'))

peaks2 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=peaksets$tvcAcc,
      denovoPeaks=setdiff(
        Reduce(union,peaksets[c('timeDep','mespDep','handrDep')]),
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)
peaks2 <- mapply(
  intersect,
  peaks2,
  list(
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.70|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .70
    ],
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.70|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.70
    ]
  )
)
sapply(peaks2,length)
peakHyper(
  peaks2,
  'mek0.7dn0.7',
  motifMatches(cisbp.matches)[
    ,names(cisbp.motifs)[
      !grepl('^pTH[0-9]+',names(cisbp.motifs))&
        !grepl("^M[A][0-9]+",names(cisbp.motifs))&
        !grepl("^Zen_Cell",names(cisbp.motifs))&
        !grepl("^ftz",names(cisbp.motifs))&
        !grepl("^Hnf4",names(cisbp.motifs))&
        !grepl("^MA[0-9]+",names(cisbp.motifs))&
        !grepl("^GM[0-9]+",names(cisbp.motifs))&
        cisbp.family!='promoter'
    ]
  ],
  cisbp.family[
      !grepl('^pTH[0-9]+',names(cisbp.motifs))&
        !grepl("^M[A][0-9]+",names(cisbp.motifs))&
        !grepl("^Zen_Cell",names(cisbp.motifs))&
        !grepl("^ftz",names(cisbp.motifs))&
        !grepl("^Hnf4",names(cisbp.motifs))&
        !grepl("^MA[0-9]+",names(cisbp.motifs))&
        !grepl("^GM[0-9]+",names(cisbp.motifs))&
        cisbp.family!='promoter'
  ],
  p=.005,fdr=F,logOR = 2.00,maskOR=T
)

peaks <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(
        Reduce(union,peaksets[c("open6",'closed6','closed18')]),
        peaksets$open18
      ),
      denovoPeaks=setdiff(
        peaksets$open18,
        Reduce(union,peaksets[c("open6",'closed6','closed18')])
      )
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks <- mapply(intersect,peaks,append(
  replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
    mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange<0
  ],F),
  replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
    mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange<0
  ],F)
))

peakHyper(
  peaks,
  "denovo.cardiac.asmFDR",
  motifMatches(matches)[
    # Reduce(union,peaksets[c("timeDep","mespDep","handrDep")]),
    unique(geneToPeak(
      con,
      union(scrna$denovoCardiac,scrna$denovoASM)
      # peaksets$timeDep
    )$PeakID),
  ],
  tf.family,p=.1,fdr = T
)


peakFamilyHyper(
  peaks,
  "denovo.cardiac.asmFamily",
  homer.family,
  motifMatches(homer.matches)[
    # Reduce(union,peaksets[c("timeDep","mespDep","handrDep")]),
    unique(geneToPeak(
      con,
      union(scrna$denovoCardiac,scrna$denovoASM)
      # peaksets$timeDep
    )$PeakID),
    fdr=F,p=.1
  ]
)

peaks2 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=peaksets$tvcAcc,
      denovoPeaks=setdiff(
        Reduce(union,peaksets[c("timeDep",'mespDep','handrDep')]),
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks2 <- mapply(intersect,peaks2,append(
  replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
    mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.8|
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .8
  ],F),
  replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
    mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.8|
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.8
  ],F)
))

peakHyper(
  peaks2,
  "denovo.cardiac.asmTVCpeakome",
  motifMatches(matches),
  tf.family,p=.05,fdr = F
)


peakHyper(
  peaks,
  "denovo.cardiac.asmVsPeakome",
  motifMatches(homer.matches),
  homer.family,
  p=.01
)


peaks2 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=peaksets$tvcAcc,
      denovoPeaks=setdiff(
        Reduce(union,peaksets[c('timeDep','mespDep','handrDep')]),
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks2 <- mapply(
  intersect,
  peaks2,
  list(
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.50
        # mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .50
    ],
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.50
        # mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.50
    ]
  )
)
sapply(peaks2,length)

peakHyper(
  peaks2,
  "denovo.cardiac.asmTVC",
  motifMatches(matches)[
    # Reduce(union,peaksets[c("timeDep","mespDep","handrDep")]),
    unique(geneToPeak(
      con,
      union(scrna$denovoCardiac,scrna$denovoASM)
      # peaksets$timeDep
    )$PeakID),
  ],
  tf.family,p=.05,fdr = F
)

peakHyper(
  peaks2,
  "denovo.cardiac.asmTVCpeakome",
  motifMatches(matches),
  tf.family,p=.05,fdr = F
)

peaks3 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoPeaks=setdiff(
        peaksets$open18,
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks3 <- mapply(
  intersect,
  peaks3,
  list(
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< 0|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> 0
    ],
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< 0|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>0
    ]
  )
)
sapply(peaks3,length)

peakHyper(
  peaks3,
  "denovo.cardiac.asmTVC_18",
  motifMatches(matches),
  tf.family,p=.05,fdr = F
)
peaks3 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=peaksets$tvcAcc,
      denovoPeaks=setdiff(
        union(peaksets$mespDep,peaksets$handrDep),
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks3 <- mapply(
  intersect,
  peaks3,
  append(
    replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< 0&
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> 0
    ],F),
    replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< 0&
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>0
    ],F)
  )
)

sapply(peaks2,length)

peakHyper(
  peaks2,
  "denovo.cardiac.asmTVCpeakome3",
  motifMatches(matches),
  tf.family,p=.05,fdr = F
)

cisbp.motifs <- getCisbpMotifs()
names(cisbp.motifs) <- make.names(sapply(tags(cisbp.motifs),'[[',"DBID.1"),T)
cisbp.matches <- matchMotifs(cisbp.motifs,ann$peaks,BSgenome.Cintestinalis.KH.KH2013,'subject','matches')
cisbp.family <- sapply(tags(cisbp.motifs),'[[',"Family_Name")

denovoMotifs <- make.names(c(
  'EBF',"V$CREB_01","JDP2_2","Jundm2_0911","Otx2_3441","RARG_3",
  "Crx_3485","Pitx2_2274","Gsc_2327","FOXF2_f1","Hoxb6_3428","TATA-box","V$OCT1_04",
  "TBX1_3","EHF_1","RORA_1",'K562_GATA2_UChicago',"SOX8_6",
  "ISL1_f1",
  "pTH4582","AP2B_f1","tgo_ss_SANGER_5_FBgn0003513","BREd",#"HeLa-S3_JUND_Stanford",
  'pTH5002',"SUH_f1",#sonic hedgehog
  "DCEI-DCEIII1",
  "Mv63",
  "MA0498.1",#MEIS1
  "Ara_Cell_FBgn0015904",
  "Mv90",#MEF2
  'pTH9356',#SAND
  "V$CDXA_02","RAX_1","V$NKX3A_01","PDX1_2","PAX6_f1",
  "pTH5709",#nuclear receptor
  "NR4A2_si",
  "RARA_6","SMAD3_f1","Med_FlyReg_FBgn0011655"
))
motifsel <- sapply(tags(cisbp.motifs),'[[',"DBID.1")%in%denovoMotifs

peakHyper(
  peaks2,
  'cisbpDenovo',
  motifMatches(cisbp.matches)[,denovoMotifs],
  cisbp.family[denovoMotifs],
  p=.05,fdr=F
)

motifsel <- cisbp.family%in%cisbp.family[denovoMotifs]

peakHyper(
  peaks,
  'cisbpDenovo18',
  motifMatches(cisbp.matches)[,denovoMotifs],
  cisbp.family[denovoMotifs],
  p=.05,fdr=F
)

peakHyper(
  peaks,
  'cisbp18',
  motifMatches(cisbp.matches),
  cisbp.family,
  p=.001,fdr=F
)

peakHyper(
  peaks3,
  'cisbpDenovoTVC18',
  motifMatches(cisbp.matches)[,denovoMotifs],
  cisbp.family[denovoMotifs],
  p=.05,fdr=F
)

peakHyper(
  peaks3,
  'cisbpTVC18',
  motifMatches(cisbp.matches),
  cisbp.family,
  p=.010,fdr=T,logOR = 1.5,maskOR = F
)

comb.motifs <- append(
  setNames(known.motifs,paste0('HOMER.',names(known.motifs))),
  cisbp.motifs[-10]
)
comb.family <- sapply(tags(comb.motifs),'[[',"Family_Name")

comb.matches <- matchMotifs(comb.motifs,ann$peaks,BSgenome.Cintestinalis.KH.KH2013,'subject','matches')

peakHyper(
  peaks3,
  'combTVC18',
  motifMatches(comb.matches),
  comb.family,
  p=.010,fdr=T,logOR = 1.5,maskOR = F
)

peaks4 <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(peaksets$tvcAcc,peaksets$open18),
      denovoPeaks=setdiff(
        peaksets$open18,
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks4 <- mapply(
  intersect,
  peaks4,
  append(
    replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< 0|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> 0
    ],F),
    replicate(2,row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< 0|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>0
    ],F)
  )
)

peakHyper(
  peaks3,
  'cisbpTVC18new',
  motifMatches(cisbp.matches),
  cisbp.family,
  p=.010,fdr=T,logOR = 1.5,maskOR = F
)


peakHyper(
  peaks2,
  'cisbpTVC18old',
  motifMatches(cisbp.matches),
  cisbp.family,
  p=.010,fdr=T,logOR = 1.5,maskOR = F
)

peaks <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=peaksets$tvcAcc,
      denovoPeaks=setdiff(
        Reduce(union,peaksets[c('timeDep','mespDep','handrDep')]),
        peaksets$tvcAcc)
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)
peaks <- mapply(
  intersect,
  peaks,
  list(
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_MekMut_vs_control)[
      mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange< -.50|
        mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange> .50
    ],
    names(ann$peaks),
    row.names(mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control)[
      mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange< -.50|
        mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange>.50
    ]
  )
)
peakFamilyHyper(
  peaks,
  'mek0.5dn0.5FamilyCombn',
  comb.family,
  motifMatches(comb.matches),
  p=.050,fdr=T,logOR = 0.00,maskOR=F
)

sapply(peaks,length)
peakHyper(
  peaks,
  'mek0.7dn0.7',
  motifMatches(cisbp.matches)[
    ,names(cisbp.motifs)[
      !grepl('^pTH[0-9]+',names(cisbp.motifs))&
        !grepl("^Zen_Cell",names(cisbp.motifs))&
        !grepl("^ftz",names(cisbp.motifs))&
        !grepl("^Hnf4",names(cisbp.motifs))&
        !grepl("^MA[0-9]+",names(cisbp.motifs))&
        !grepl("^GM[0-9]+",names(cisbp.motifs))&
        cisbp.family!='promoter'
    ]
  ],
  cisbp.family[
      !grepl('^pTH[0-9]+',names(cisbp.motifs))&
        !grepl("^Zen_Cell",names(cisbp.motifs))&
        !grepl("^ftz",names(cisbp.motifs))&
        !grepl("^Hnf4",names(cisbp.motifs))&
        !grepl("^MA[0-9]+",names(cisbp.motifs))&
        !grepl("^GM[0-9]+",names(cisbp.motifs))&
        cisbp.family!='promoter'
  ],
  p=.010,fdr=T,logOR = 2.00,maskOR=F
)

peakFamilyHyper(
  peaks,
  'mek0.7dn0.7Family',
  cisbp.family,
  motifMatches(cisbp.matches),
  p=.050,fdr=T,logOR = 0.00,maskOR=F
)

peakFamilyHyper(
  peaks,
  'mek0.7dn0.7FamilyHomer',
  tf.family,
  motifMatches(matches),
  p=.050,fdr=T,logOR = 0.00,maskOR=F
)


sel <- !grepl('^pTH[0-9]+',names(cisbp.motifs))&
        !grepl("^Zen_Cell",names(cisbp.motifs))&
        !grepl("^Caup_Cell",names(cisbp.motifs))&
        !grepl("^ftz",names(cisbp.motifs))&
        !grepl("^Hnf4",names(cisbp.motifs))&
        !grepl("^M[Afvw][0-9]+",names(cisbp.motifs))&
        !grepl("^GM[0-9]+",names(cisbp.motifs))&
        cisbp.family!='promoter'
peakHyper(
  peaks,
  'TVC_18',
  motifMatches(cisbp.matches)[
    ,names(cisbp.motifs)[sel]
  ],
  cisbp.family[sel],
  p=.050,fdr=T,logOR = 1.70,maskOR=F,breaks = c(1,4)
)

peakFamilyHyper(
  peaks,
  'homer_TVC_18Family',
  tf.family,
  motifMatches(matches),
  p=.050,fdr=T,logOR = 0.00,maskOR=T
)
