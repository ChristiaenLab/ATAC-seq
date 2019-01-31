#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=30GB
#SBATCH --output=counts_%a.out
#SBATCH --error=counts_%a.err

module load bedtools/intel/2.26.0

cd /scratch/kaw504/atacCiona/atac/bam
out="/scratch/kaw504/atacCiona/DE_DA"

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" file_name.txt )
mkdir -p ${out}/counts.new

bedtools multicov -bams  ${sample}_q30_rmdup_KhM0_sorted.bam -bed ${out}/peakome.bed > ${out}/counts.new/${sample}_counts_new.txt

exit 0;
