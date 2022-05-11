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

for j in {1..10}; do
python ~/software/PRScsx/PRScsx.py \
  --ref_dir=../../../../ref/PRS_CSx/ \
  --bim_prefix=../../../../genotype/EUR/test/test_5k \
  --sst_file=../../../../summary_stat/${pops}/Scene1_20K/${pop1}/discover/sim_${i}_PRScs.txt,../../../../summary_stat/${pops}/Scene1_20K/${pop2}/discover/sim_${i}_PRScs.txt \
  --n_gwas=40000,20000 \
  --pop=${pop1},${pop2} \
  --out_dir=./ \
  --out_name=sim_${i} \
  --chrom=${j}

module load PLINK/1.90-beta5.3

mkdir -p ../predict/PRS_CSx

#plink --bfile ../../../../genotype/EUR/test/test_5k --score sim_${i}_EUR_pst_eff_a1_b0.5_phiauto_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EUR_auto
i
#plink --bfile ../../../../genotype/EAS/test/test_5k --score sim_${i}_EAS_pst_eff_a1_b0.5_phiauto_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EAS_auto

for phi in 1e-6 1e-4 1e-2 1;
do
    mkdir -p ${phi}/
    python ~/software/PRScsx/PRScsx.py --ref_dir=../../../../ref/PRS_CSx/ --bim_prefix=../../../../genotype/EUR/test/test_5k --sst_file=../../../../summary_stat/${pops}/Scene1_20K/${pop1}/discover/sim_${i}_PRScs.txt,../../../../summary_stat/${pops}/Scene1_20K/${pop2}/discover/sim_${i}_PRScs.txt --n_gwas=40000,20000 --pop=${pop1},${pop2} --out_dir=${phi}/ --out_name=sim_${i} --phi=${phi} --chrom=${j}

    #mv ${phi}/sim_${i}_EUR_*.txt ${phi}/sim_${i}_EUR_chr${j}.txt
    #mv ${phi}/sim_${i}_EAS_*.txt ${phi}/sim_${i}_EAS_chr${j}.txt

    #plink --bfile ../../../../genotype/EUR/test/test_5k --score ${phi}/sim_${i}_EUR_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EUR_${phi}
    #plink --bfile ../../../../genotype/EAS/test/test_5k --score ${phi}/sim_${i}_EAS_chr1.txt 2 4 6 --out ../predict/PRS_CSx/sim_${i}_EAS_${phi}
done
done

rm -f ../predict/PRS_CSx/*.nosex
