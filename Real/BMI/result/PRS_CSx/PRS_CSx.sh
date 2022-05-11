#!/bin/bash
#SBATCH -J PRS_CSx
#SBATCH --partition scavenge 
#SBATCH -c 3
#SBATCH --mem 10g
#SBATCH --time=1-00:00:00 

##################################### EUR_EAS ##################################
pop1=EUR
pop2=EAS
pops=EUR_EAS

i=${SLURM_ARRAY_TASK_ID}

source activate python2

python ../../../software/PRScsx/PRScsx.py \
  --ref_dir=/ysm-gpfs/pi/zhao/gz222/UKB_real/ref/PRS_CS/ \
  --bim_prefix=../../../genotype/EAS/Ukb_imp_v2 \
  --sst_file=/ysm-gpfs/pi/zhao-data/gz222/height/summ_stats/PRS_cs.txt,../../summary_stat/EAS/PRS_cs.txt \
  --n_gwas=233766,158284 \
  --pop=${pop1},${pop2} \
  --out_dir=./ \
  --out_name=PRS_cs \
  --chrom=${i}

module load PLINK/1.90-beta5.3

mkdir -p ../predict/PRS_CSx

#plink --bfile ../../../../genotype/EUR/test/test_5k --score sim_${i}_EUR_pst_eff_a1_b0.5_phiauto_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EUR_auto
#plink --bfile ../../../../genotype/EAS/test/test_5k --score sim_${i}_EAS_pst_eff_a1_b0.5_phiauto_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EAS_auto

for phi in 1e-6 1e-4 1e-2 1;
do
    mkdir -p ${phi}/
    python ../../../software/PRScsx/PRScsx.py --ref_dir=/ysm-gpfs/pi/zhao/gz222/UKB_real/ref/PRS_CS/ --bim_prefix=../../../genotype/EAS/Ukb_imp_v2 --sst_file=/ysm-gpfs/pi/zhao-data/gz222/height/summ_stats/PRS_cs.txt,../../summary_stat/EAS/PRS_cs.txt --n_gwas=233766,158284 --pop=${pop1},${pop2} --out_dir=${phi}/ --out_name=PRS_cs --phi=${phi} --chrom=${i}

    #mv ${phi}/sim_${i}_EUR_*.txt ${phi}/sim_${i}_EUR_chr${j}.txt
    #mv ${phi}/sim_${i}_EAS_*.txt ${phi}/sim_${i}_EAS_chr${j}.txt

    #plink --bfile ../../../../genotype/EUR/test/test_5k --score ${phi}/sim_${i}_EUR_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EUR_${phi}
    #plink --bfile ../../../../genotype/EAS/test/test_5k --score ${phi}/sim_${i}_EAS_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EAS_${phi}
done

rm -f ../predict/PRS_CSx/*.nosex
