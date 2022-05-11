#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 15g
#SBATCH --time 05:00:00
#SBATCH --job-name ldpred2_ref

i=${SLURM_ARRAY_TASK_ID}

module load R

Rscript calc_ref.R ${i} 
