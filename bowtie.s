#!/bin/bash
#
#SBATCH --job-name=bowtie
#SBATCH --nodes=1 --tasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mem=62GB
#SBATCH --output=bowtie%a.out
#SBATCH --error=bowtie%a.err

module purge
module load bowtie2/intel/2.3.2
module load samtools/intel/1.3.1
module load ngsutils/intel/0.5.9
module load picard/2.8.2
module load preseq/intel/2.0.1
module load deeptools/intel/2.4.2
module load  bedtools/intel/2.26.0
module load r/intel/3.4.2
module load trim_galore/0.4.4
module load cutadapt/intel/1.12


RUNDIR=/scratch/kaw504/atacCiona/atac

#cd /scratch/kaw504/atacCiona/august
#bowtie2-build -f JoinedScaffold JoinedScaffold

cd $RUNDIR
N1=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" bowtie/fastq.txt)

N2=${N1/n01/n02}
sample=${N1/n01/}

#bowtie2 --no-discordant -p 20 --no-mixed -N 1 -X 1000  --un-conc bowtie/${sample}.unaligned.fastq -x /scratch/kaw504/atacCiona/august/JoinedScaffold -1 fastq_files/atac/${N1}*.fastq.gz  -2 fastq_files/atac/${N2}*.fastq.gz -S bowtie/${sample}.sam 

cd bowtie

#samtools view -bS -q 30 ${sample}.sam > ${sample}_q30.bam
#rm ${sample}.sam

trim_galore --nextera --paired --gzip --three_prime_clip_R1 1 --three_prime_clip_R2 1 ${sample}.unaligned.1.fastq ${sample}.unaligned.2.fastq

# Align trimmed reads
bowtie2 --no-discordant -p 20 --no-mixed -N 1 -X 1000  --un-conc ${sample}.still_unaligned.fastq -x /scratch/kaw504/atacCiona/august/JoinedScaffold -1 ${sample}.unaligned.1_val_1.fq.gz  -2 ${sample}.unaligned.2_val_2.fq.gz -S ${sample}.aligned.sam
# Convert trimmed read results to BAM
samtools view -Sb -q 30 ${sample}.aligned.sam > ${sample}.trimmed.bam
rm ${sample}.aligned.sam

# Concatenate initial and trimmed results
samtools cat -o ${sample}.merged.bam ${sample}_q30.bam ${sample}.trimmed.bam
# Retain only aligned pairs based on BAM flags
samtools view -b -F 4 ${sample}.merged.bam > ${sample}.mapped.bam
# Sort mapped reads
samtools sort -o ${sample}.srt.bam ${sample}.mapped.bam 
# Remove duplicates
samtools rmdup ${sample}.srt.bam ${sample}.rmd.srt.bam 

java -Xmx10g -jar /share/apps/picard/2.8.2/picard-2.8.2.jar MarkDuplicates INPUT=${sample}.rmd.srt.bam OUTPUT=${sample}_q30_rmdup_sorted.bam METRICS_FILE=${sample}_DEDUPL.out.mat REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2048

bamutils filter ${sample}_q30_rmdup_sorted.bam  ${sample}_q30_rmdup_KhM0_sorted.bam -excluderef KhM0
# Index the unique, sorted reads
samtools index ${sample}_q30_rmdup_KhM0_sorted.bam
# Compute stats for all mapped reads and for reads after removing duplicates
samtools flagstat ${sample}.merged.bam > ${sample}.stats
samtools flagstat ${sample}.rmd.srt.bam >> ${sample}.stats
samtools flagstat ${sample}_q30_rmdup_KhM0_sorted.bam >> ${sample}.stats

# Perform preseq on sorted reads to assess library diversity
preseq lc_extrap -P -o ${sample}.txt <( bamToBed -i ${sample}.srt.bam )
# Perform insert size analysis on unique, sorted reads
java -Xmx10g -jar /share/apps/picard/2.8.2/picard-2.8.2.jar CollectInsertSizeMetrics I=${sample}_q30_rmdup_KhM0_sorted.bam O=${sample}_insert_sizes.txt H=${sample}.pdf

bamCoverage -of bigwig -o ${sample}.bw -b ${sample}_q30_rmdup_KhM0_sorted.bam --ignoreDuplicates -p 20 --minMappingQuality 30 --centerReads -ignore KhM0 --normalizeUsingRPKM

#igvtools count -w 50 ${sample}_q30_rmdup_KhM0_sorted.bam ${sample}_q30_rmdup_KhM0_sorted.tdf /scratch/cr1636/august/JoinedScaffold 

#bam2wig.py --input-file=${sample}_q30_rmdup_KhM0_sorted.bam --chromSize=/scratch/kaw504/atacCiona/august/JoinedScaffold.fai --out-prefix=${sample}_q30_rmdup_KhM0_sorted

#wigToBigWig ${sample}_q30_rmdup_KhM0_sorted.wig /scratch/kaw504/atacCiona/august/JoinedScaffold.fai ${sample}_q30_rmdup_KhM0_sorted.bw


exit 0;

