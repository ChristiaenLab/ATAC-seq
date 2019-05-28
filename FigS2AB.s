#!/bin/bash
#
#SBATCH --job-name=plotprofile
#SBATCH --nodes=1 --tasks-per-node=20
#SBATCH --time=1:00:00
#SBATCH --mem=62GB
#SBATCH --output=plotprofile.out
#SBATCH --error=plotprofile.err


module purge
module load bedtools/intel/2.26.0
module load deeptools/intel/2.4.2
RUNDIR=/scratch/kaw504/atacCiona/DE_DA

cd $RUNDIR

bedtools intersect -f 1 -a foxf.cisbp.bed -b 2018-10-03/tvcAcc.bed > foxf.tvc.bed

computeMatrix reference-point -S \
 /scratch/kaw504/atacCiona/atac/bam/3xbfp_5hpf_26_july.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_lacz_10hpf_21_april.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_mekmut_10hpf_14_june.bw \
 /scratch/kaw504/atacCiona/atac/bam/foxf_22_2.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_dnfgfr_10hpf_26_april.bw \
 -R foxf.tvc.bed --skipZeros -p 20 -o foxf.tvc.mat.gz --referencePoint center -b 500 -a 500

plotProfile -m foxf.tvc.mat.gz --outFileName foxf.tvc.eps --perGroup --legendLocation lower-right --refPointLabel "FoxF binding site" --samplesLabel control "mesp MekMut" "FoxF KO"

computeMatrix reference-point -S \
 /scratch/kaw504/atacCiona/atac/bam/3xbfp_5hpf_26_july.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_lacz_10hpf_21_april.bw \
 /scratch/kaw504/atacCiona/atac/bam/lacz1_22_3.bw \
 /scratch/kaw504/atacCiona/atac/bam/lacz_8_3.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_lacz_10hpf_2_june.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_dnfgfr_10hpf_26_april.bw \
 /scratch/kaw504/atacCiona/atac/bam/mesp_mekmut_10hpf_14_june.bw \
 /scratch/kaw504/atacCiona/atac/bam/foxf_22_2.bw \
 -R foxf.cisbp.bed --skipZeros -p 20 -o foxf.cisbp.mat.gz --referencePoint center -b 500 -a 500
# /scratch/kaw504/atacCiona/atac/bam/3xbfp_5hpf_16_aug.bw \
# /scratch/kaw504/atacCiona/atac/bam/3xbfp_5hpf_21_july.bw \
# /scratch/kaw504/atacCiona/atac/bam/mesp_dnfgfr_10hpf_23_3_17.bw \
# /scratch/kaw504/atacCiona/atac/bam/mesp_mekmut_10hpf_3_june_2016.bw \
# /scratch/kaw504/atacCiona/atac/bam/foxf_23_3.bw \

plotProfile -m foxf.cisbp.mat.gz --outFileName foxf.eps --perGroup --legendLocation lower-right --refPointLabel "FoxF binding site" #--samplesLabel control "mesp dnFGFR" "mesp MekMut" "FoxF KO"
exit 0;

