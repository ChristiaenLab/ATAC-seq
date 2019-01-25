#!/bin/bash

# submits an R script as a batch job
# usage: ./rscript.s <script.R> <optional arguments>
echo "#!/bin/bash
#
#SBATCH --job-name=rscript
#SBATCH --cpus-per-task=20
#SBATCH --time=3:00:00
#SBATCH --mem=62GB
#SBATCH --output=r%A.out
#SBATCH --error=r%A.err 

module purge
module load gcc/6.3.0
module load gsl/intel/2.3
module load r/intel/3.4.2

Rscript $@

exit 0;
" | sbatch

