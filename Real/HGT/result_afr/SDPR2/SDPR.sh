#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 15g
#SBATCH --time 2-00:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}

source activate python2

python2 ../../../software/SDPRX2/SDPRX.py --ss1 /ysm-gpfs/pi/zhao-data/gz222/height/summ_stats/SDPR.txt --ss2 ../../summary_stat/AFR/SDPR.txt --N1 252230 --N2 49781 --load_ld ../../../ref/AFR/SDPRX/ --valid ../../../genotype/AFR/Ukb_imp_v2.bim --chr ${i} --rho 0.86 --out SDPR_${i}

#cat *.out | grep -oP '(?<=h2: ).*' > h2.txt

#module load PLINK/1.90-beta5.3
#plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score res.txt 1 2 3 --out ../predict/SDPR/res
