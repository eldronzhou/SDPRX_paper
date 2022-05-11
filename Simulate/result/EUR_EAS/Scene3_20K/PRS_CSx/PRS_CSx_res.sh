#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 3g
#SBATCH --time 1:00:00
#SBATCH --job-name PRS_cSx

i=${SLURM_ARRAY_TASK_ID}
j=Scene3
k=50K

cat sim_${i}_EUR*.txt > sim_${i}_EUR.txt
cat 1e-6/sim_${i}_EUR*.txt > 1e-6/sim_${i}_EUR.txt
cat 1e-4/sim_${i}_EUR*.txt > 1e-4/sim_${i}_EUR.txt
cat 1e-2/sim_${i}_EUR*.txt > 1e-2/sim_${i}_EUR.txt
cat 1/sim_${i}_EUR*.txt > 1/sim_${i}_EUR.txt

cat sim_${i}_EAS*.txt > sim_${i}_EAS.txt
cat 1e-6/sim_${i}_EAS*.txt > 1e-6/sim_${i}_EAS.txt
cat 1e-4/sim_${i}_EAS*.txt > 1e-4/sim_${i}_EAS.txt
cat 1e-2/sim_${i}_EAS*.txt > 1e-2/sim_${i}_EAS.txt
cat 1/sim_${i}_EAS*.txt > 1/sim_${i}_EAS.txt

module load R
mkdir -p ../predict/PRS_CSx
h2=0.3
Rscript cv.R ${i} ${h2} ${j}

rm -f sim_${i}_*
rm -f 1e-6/sim_${i}_*
rm -f 1e-4/sim_${i}_*
rm -f 1e-2/sim_${i}_*
rm -f 1/sim_${i}_*



