#!/bin/bash
#
#SBATCH --job-name=star
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH --output=star%a.out
#SBATCH --error=star%a.err

module purge
module load star/intel
module load samtools/intel/1.3.1
module load ngsutils/intel
module load htseq/intel

RUNDIR=/scratch/cr1636/ATAC_Ciona_sana_claudia/RNAseq_foxf_ngn_lacz_june2017/2017-06-28_HKHMWAFXX/merged/

cd $RUNDIR

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" lib.txt)

STAR --runThreadN 12 --genomeDir star --readFilesIn HKHMWAFXX_n01_${sample}.fastq HKHMWAFXX_n02_${sample}.fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts TranscriptomeSAM --bamRemoveDuplicatesType UniqueIdentical --outWigReferencesPrefix ${sample} --outWigType wiggle --outFileNamePrefix ${sample}

RUNDIR=/scratch/cr1636/ATAC_Ciona_sana_claudia/RNAseq_foxf_ngn_lacz_june2017/2017-06-28_HKHMWAFXX/merged/

cd $RUNDIR

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" lib.txt)

samtools rmdup ${sample}Aligned.sortedByCoord.out.bam ${sample}Aligned.sortedByCoord.out_rmdup.bam 
bamutils filter ${sample}Aligned.sortedByCoord.out_rmdup.bam  ${sample}_dupAligned.sortedByCoord.out_rmdup_KhM0.bam -excluderef KhM0
samtools index -b ${sample}Aligned.sortedByCoord.out_rmdup_KhM0.bam

htseq-count -f bam -r pos -i gene_id ${sample}Aligned.sortedByCoord.out_rmdup_KhM0.bam  /scratch/cr1636/august/KH.KHGene.2013.gtf > ${sample}_star_count.txt

exit 0;

