#!/bin/bash
#
#SBATCH --job-name=macs2_broad
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --output=wordcounts_%A_%a.out
#SBATCH --error=wordcounts_%A_%a.err


module load macs2/intel/2.1.1

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" file_name.txt)

macs2 callpeak --nomodel -t ${sample}_q30_rmdup_KhM0_sorted.bam -c gdna_4jun_q30_rmdup_KhM0_sorted.bam -f BAM -g 1.15e8 --keep-dup all --shift -29 --extsize 50 --outdir ./macs2_shift29_ext50_withINPUT  -n ${sample} --bdg --call-summits 

exit 0;

