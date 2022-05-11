#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 15g
#SBATCH --time 1-00:00:00
#SBATCH --job-name ldpred2

i=${SLURM_ARRAY_TASK_ID}
j=Scene2

module load R

mkdir -p sim_${i}/EAS

export OPENBLAS_NUM_THREADS=1

Rscript ldpred2_pipeline.R ${i} ${j}_20K

rm -f sim_${i}/EAS/.sbk

mkdir -p ../predict/ldpred2

Rscript cv.R ${i} 0.3 ${j}

rm -rf sim_${i}
