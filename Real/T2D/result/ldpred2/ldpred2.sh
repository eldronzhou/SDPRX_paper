#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 10
#SBATCH --mem 35g
#SBATCH --time 1-00:00:00
#SBATCH --job-name ldpred2

module load R

export OPENBLAS_NUM_THREADS=1

Rscript ldpred2_pipeline.R 

mkdir -p ../predict/ldpred2

module load PLINK

plink2 --bfile ../../../genotype/EAS/Ukb_imp_v2 --score tmp.txt 1 2 --score-col-nums 3-67 --threads 3 --memory 5120 --out ../predict/ldpred2/ldpred2

rm -rf tmp-data/
gzip ../predict/ldpred2/ldpred2.sscore
