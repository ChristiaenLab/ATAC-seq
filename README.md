# Ciona ATAC-seq data analysis pipeline
This pipeline is designed for processing of paired-end ATAC-seq libraries.
The pipeline can be run on compute clusters with job submission engines or stand alone machines. Beginning from raw FASTQ files, the pipeline calls peaks and generates signal tracks. An accessome is created from all samples, which is used to compute read counts to calculate differential accessibility. The pipeline integrates ATAC-seq and RNA-seq data by annotating peaks to nearby genomic elements, and merging differentially expressed genes with differentially accessible peaks. It further characterizes differentially accessible elements by performing Gene Set Enrichment Analysis and motif enrichment.

The shell scripts (*.s) are intended for submission using slurm, and may require modification before they can be run with other job submission managers.
rscript.s is a wrapper script for submitting R scripts as batch jobs. All R scripts (*.R) can instead be run interactively with Rscript.

Tools required: bamutils, bedtools, Bowtie2, deeptools, gcc>=6.3, gsl>=2.3, htseq, kent, macs2, picard, R>=3.4 

R packages required:  BSgenomes.Cintestinalis.KH.JoinedScaffold, circlize, chromVAR, ComplexHeatmap, DBI, DESeq2, edgeR, fgsea, GenomicFeatures, GenomicRanges, ggplot2, lattice, latticeExtra, motifmatchr, optparse, reshape2, RSQLite, rtracklayer, TFBSTools, UpSetR, VennDiagram

----------------------------
Genome data

Usage:

get KH2013 transcript model from GHOST
```bash
wget -U firefox http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2013.gff3.zip
unzip KH.KHGene.2013.gff3.zip
```

install BSgenomes.Cintestinalis.KH.JoinedScaffold
```bash
cd BSgenome.Cintestinalis.KH.JoinedScaffold
./forgeGenome.sh
```
align ATACseq reads
```bash
sbatch -a1-52 bowtie.s
```
perform QC, generate bigWig files for IGV plots
```bash
sbatch -a1-52 seqstats.s
```
call peaks
```bash
sbatch -a1-52 macs2.s
```
merge peaks into accessome
```bash
sbatch peakome.s
```
get read counts for peaks in accessome
```bash
sbatch -a1-52 counts.s
```
create STAR index
```bash
sbatch staridx.s
```
align RNAseq reads and get counts
```bash
sbatch -a1-6 star.s
```
annotate peaks and initialize database
```bash
./rscript.s writeDB.R
```
calculate differential expression & add to database
```bash
./rscript.s rnaseqFoxf.R
./rscript.s rnaseqMAPK.R
```
calculate differential accessibility & add to database
```bash
./rscript.s atacDESeq.R
```
perform GSEA on peak accessibility
```bash
./rscript.s runFgsea.R
```
motif analysis of DA peaks
```bash
./rscript runChromVar.R
```
