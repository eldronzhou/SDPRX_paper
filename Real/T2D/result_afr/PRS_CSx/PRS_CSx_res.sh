#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 1-00:00:00
#SBATCH --job-name PRS_cS

cat PRS_cs_EUR_*.txt > PRS_cs_EUR.txt
cat PRS_cs_AFR_*.txt > PRS_cs_AFR.txt

for phi in 1e-6 1e-4 1e-2 1
do
cat ${phi}/PRS_cs_EUR_*.txt >${phi}/PRS_cs_EUR.txt
cat ${phi}/PRS_cs_AFR_*.txt > ${phi}/PRS_cs_AFR.txt
done

module load PLINK/1.90-beta5.3

mkdir -p ../predict/PRS_CS/

plink --bfile ../../../genotype/AFR/Ukb_imp_v2 --score PRS_cs_AFR.txt 2 4 6 --out ../predict/PRS_CS/PRS_CS_AFR

for phi in 1e-6 1e-4 1e-2 1
do
    plink --bfile ../../../genotype/AFR/Ukb_imp_v2 --score ${phi}/PRS_cs_AFR.txt 2 4 6 --out ../predict/PRS_CS/PRS_CS_AFR_${phi}
done

gzip ../predict/PRS_CS/*.profile


