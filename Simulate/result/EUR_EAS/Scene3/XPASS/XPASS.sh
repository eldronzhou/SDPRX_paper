#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 10g
#SBATCH --time 1:00:00
#SBATCH --job-name XPASS

i=${SLURM_ARRAY_TASK_ID}
k=Scene3

module load R
module load PLINK/1.90-beta5.3

Rscript XPASS.R ${i} ${k}

rm -f _PosteriorMean.txt 

mkdir -p ../predict/XPASS

plink --bfile ../../../../genotype/EUR/test/test_5k --score sim_${i}.txt 1 2 3 header --out ../predict/XPASS/sim_${i}_EUR

plink --bfile ../../../../genotype/EAS/test/test_5k --score sim_${i}.txt 1 2 4 header --out ../predict/XPASS/sim_${i}_EAS
