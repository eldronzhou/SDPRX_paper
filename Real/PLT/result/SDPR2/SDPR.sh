#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 15g
#SBATCH --time 1-00:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}

source activate python2

python2 ../../../software/SDPRX2/SDPRX.py --ss1 ../../summary_stat/EUR/SDPR.txt --ss2 ../../summary_stat/EAS/SDPR.txt --N1 563085 --N2 108208 --load_ld ../../../ref/EAS/SDPRX/ --valid ../../../genotype/EAS/Ukb_imp_v2.bim --chr ${i} --threads 3 --rho 0.84 --out SDPR_${i}

#cat *.out | grep -oP '(?<=h2: ).*' > h2.txt

#module load PLINK/1.90-beta5.3
#plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score res.txt 1 2 3 --out ../predict/SDPR/res
