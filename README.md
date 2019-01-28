# Ciona ATAC-seq data analysis pipeline
This pipeline is designed for processing of paired-end ATAC-seq libraries.
It can be run on compute clusters with job submission engines or stand alone machines. The pipeline can be run starting from raw FASTQ files all the way to peak calling and signal track generation.

Tools required: bamutils, bedtools, Bowtie2, deeptools, gcc>=6.3, gsl>=2.3, htseq, macs2, R>=3.4 

R packages required:  BSgenomes.Cintestinalis.KH.KH2013, circlize, chromVAR, ComplexHeatmap, DBI, DESeq2, fgsea, GenomicFeatures, GenomicRanges, ggplot2, lattice, latticeExtra, motifmatchr, optparse, reshape2, RSQLite, rtracklayer, TFBSTools, UpSetR, VennDiagram

----------------------------
Genome data

