#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 30g
#SBATCH --time 1-00:00:00
#SBATCH --job-name XPASS

module load R
module load PLINK/1.90-beta5.3

Rscript XPASS.R 

rm -f _PosteriorMean.txt 

mkdir -p ../predict/XPASS

plink --bfile ../../../genotype/EAS/Ukb_imp_v2 --score xpass.txt header 1 2 3 --out ../predict/XPASS/EAS

plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score xpass.txt header 1 2 4 --out ../predict/XPASS/EUR
