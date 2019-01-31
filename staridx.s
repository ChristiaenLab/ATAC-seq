#!/bin/bash
#
#SBATCH --job-name=staridx
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=10:00:00
#SBATCH --mem=8GB
#SBATCH --output=staridx.out
#SBATCH --error=staridx.err

module purge
module load star/intel

RUNDIR=/scratch/cr1636/ATAC_Ciona_sana_claudia/RNAseq_foxf_ngn_lacz_june2017/2017-06-28_HKHMWAFXX/merged/

cd $RUNDIR

#sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" lib.txt)

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir star --genomeFastaFiles /scratch/cr1636/august/JoinedScaffold --sjdbGTFfile /scratch/cr1636/august/KH.KHGene.2013.gtf

exit 0;

