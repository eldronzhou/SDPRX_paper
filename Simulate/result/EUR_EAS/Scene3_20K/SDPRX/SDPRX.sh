#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 4
#SBATCH --mem 10g
#SBATCH --time 1-00:00:00
#SBATCH --job-name SDPRX

i=${SLURM_ARRAY_TASK_ID}
k=Scene3_20K

source activate python2

python ../../../../software/SDPRX.py --ss1 ../../../../summary_stat/EUR_EAS/${k}/EUR/discover/sim_${i}.PHENO2.glm.linear --ss2 ../../../../summary_stat/EUR_EAS/${k}/EAS/discover/sim_${i}.PHENO2.glm.linear --load_ld ../../../../ref/EUR_EAS/SDPRX/EUR_EAS_r.1.gz --N1 40000 --N2 20000 --mcmc_samples 1000 --rho 0 --burn 200 --VS True --out sim_${i} --threads 4

module load PLINK/1.90-beta5.3

mkdir -p ../predict/SDPRX

plink --bfile ../../../../genotype/EAS/test/test_5k --score sim_${i}_2.txt 1 2 4 header --out ../predict/SDPRX/sim_${i}_EAS

rm -f ../predict/SDPRX/*.nosex


