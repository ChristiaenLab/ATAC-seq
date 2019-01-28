library(motifmatchr)
library(TFBSTools)
library(DBI)

source('data/plotMotifs.R')

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

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

peaks <- unlist(lapply(
  scrna[c("denovoCardiac","denovoASM")],
  function(x) lapply(
    list(
      primedPeaks=setdiff(
        Reduce(union,peaksets[c("open6",'closed6','closed18')]),
        peaksets$open18
      ),
      denovoPeaks=peaksets$open18
    ),
    function(y) unique(mergeGenePeak(con,x,y)$PeakID)
  )
),F)

peaks <- mapply(intersect,peaks,append(
  replicate(2,row.names(mekmut.dnfgfr.18)[mekmut.dnfgfr.18$condition_handr_MekMut_vs_control$log2FoldChange<0],F),
  replicate(2,row.names(mekmut.dnfgfr.18)[mekmut.dnfgfr.18$condition_handr_dnFGFR_vs_control$log2FoldChange<0],F)
))

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

peakHyper(
  peaks,
  "denovo.cardiac.asm",
  motifMatches(homer.matches)[
    # Reduce(union,peaksets[c("timeDep","mespDep","handrDep")]),
    unique(geneToPeak(
      con,
      union(scrna$denovoCardiac,scrna$denovoASM)
      # peaksets$timeDep
    )$PeakID),
  ],
  homer.family,p=.1
)

peakHyper(
  peaks,
  "denovo.cardiac.asmVsPeakome",
  motifMatches(homer.matches),
  homer.family,
  p=.01
)

