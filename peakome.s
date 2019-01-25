#!/bin/bash
#
#SBATCH --job-name=ann_10kb_peakome
#SBATCH --nodes=1 
#SBATCH --time=10:00:00
#SBATCH --mem=80GB
#SBATCH --output=peakome.out
#SBATCH --error=peakome.err

module purge
module load bedtools/intel/2.26.0

cd /scratch/kaw504/atacCiona/atac/bam/macs2_shift29_ext50_withINPUT

cat \
3xbfp_5hpf_16_aug_peaks.narrowPeak \
3xbfp_5hpf_21_july_peaks.narrowPeak \
3xbfp_5hpf_26_july_peaks.narrowPeak \
dn_15hpf_14_apr_peaks.narrowPeak \
dn_15hpf_30_apr_peaks.narrowPeak \
dn_18_4_june_peaks.narrowPeak \
dn_20hpf_9_11_june_peaks.narrowPeak \
ef1_alpha_16_june_peaks.narrowPeak \
ef1_alpha_2_june_peaks.narrowPeak \
foxf_10_3_peaks.narrowPeak \
foxf_22_2_peaks.narrowPeak \
foxf_23_3_peaks.narrowPeak \
gfp_10hpf_28_april_peaks.narrowPeak \
gfp1_18hpf_12_8_peaks.narrowPeak \
gfp_15hpf_23_04_15_peaks.narrowPeak \
gfp_15hpf_30apr_peaks.narrowPeak \
gfp18_4_june_2017_peaks.narrowPeak \
gfp2_18hpf_12_08_peaks.narrowPeak \
gfp21hpf_16_june_2015_peaks.narrowPeak \
gfp_5000_27_2_peaks.narrowPeak \
gfp_6000_22_3_peaks.narrowPeak \
handr_dnfgfr_18hpf_10_6_peaks.narrowPeak \
handr_dnfgfr_18hpf_12_9_600cells_peaks.narrowPeak \
handr_lacz_18hpf_10_6_peaks.narrowPeak \
handr_lacz_18hpf_12_08_peaks.narrowPeak \
handr_mekmut1_18hpf_12_13_peaks.narrowPeak \
handr_mekmut_18hpf_10_4_peaks.narrowPeak \
handr_mekmut2_18hpf_12_13_peaks.narrowPeak \
lacz1_22_3_peaks.narrowPeak \
lacz15hpf_14_23_april_peaks.narrowPeak \
lacz_15hpf_30apr_peaks.narrowPeak \
lacz_18hpf_14_may_peaks.narrowPeak \
lacz_18hpf_19_may_4_june_peaks.narrowPeak \
lacz21hpf_16_june_2015_peaks.narrowPeak \
lacz_8_3_peaks.narrowPeak \
mesp_dnfgfr_10hpf_22_3_17_peaks.narrowPeak \
mesp_dnfgfr_10hpf_23_3_17_peaks.narrowPeak \
mesp_dnfgfr_10hpf_26_april_peaks.narrowPeak \
mesp_lacz_10hpf_21_april_peaks.narrowPeak \
mesp_lacz_10hpf_2_june_peaks.narrowPeak \
mesp_mekmut_10hpf_14_june_peaks.narrowPeak \
mesp_mekmut_10hpf_3_june_2016_peaks.narrowPeak \
mesp_mekmut_10hpf_3_may_2016_june_aug_peaks.narrowPeak \
ngn1_23_3_peaks.narrowPeak \
ngn_22_3_peaks.narrowPeak \
ngn_9_3_peaks.narrowPeak \
> peakome/peaks.bed

sort -k1,1 -k2,2n peakome/peaks.bed > peakome/peaks.sorted.bed

bedtools merge -i peakome/peaks.sorted.bed > peakome/peaks.merged.bed
bedtools subtract -A -a peakome/peaks.merged.bed -b peakome/seq_to_remove2.bed > peakome/peaks.subtracted.bed

exit 0;
