#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 2-00:00:00
#SBATCH --job-name PRS_cS

i=${SLURM_ARRAY_TASK_ID}
j=Scene3
k=10K

module load R
module load PLINK/1.90-beta5.3

mv sim_${i}_EUR_*.txt > sim_${i}_EUR_chr${j}.txt
mv sim_${i}_EAS_*.txt > sim_${i}_EAS_chr${j}.txt

for phi in 1e-6 1e-4 1e-2 1;
do
    mv ${phi}/sim_${i}_EUR_*.txt ${phi}/sim_${i}_EUR_chr${j}.txt
    mv ${phi}/sim_${i}_EAS_*.txt ${phi}/sim_${i}_EAS_chr${j}.txt
done

mkdir -p ../predict/PRS_CSx
h2=0.3
Rscript cv.R ${i} ${h2} ${j}
