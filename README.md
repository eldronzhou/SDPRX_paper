# SDPRX_paper
We are not able to release data used in the analysis because of the restricted access to the UK Biobank genotype and phenotype data. However, if you have access to UK Biobank, you can follow the instructions to get the results of this paper. The path to the software and datasets in these scripts may not be correct. We use slurm to schedule jobs to HPC, and you need to change the header if you use another system. If you have issues about these scripts, please report to the [issue](https://github.com/eldronzhou/SDPRX_paper/issues) page.

# Table of Contents
- [Simulations](#sim)
  - [Requirements](#sim-req)
  - [Obtaining genotype data](#sim-geno)
  - [Running the analysis](#sim-analysis)
- [Real data applications](#real)
  - [Requirements](#real-req)
  - [Obtaining genotype data](#real-geno)
  - [Constructing the reference LD matrix](#real-ref)
  - [Obtaining and cleaning the summary statistics](#real-ss)
  - [Running the analysis](#real-analysis)


# <a name="sim"></a>Simulations
## <a name="sim-req"></a>Requirements
* PRS-CSx 
* LDpred2 
* PLINK-1.90 
* R
* about 35 Gb desk space 

## Workflow
**<a name="sim-geno"></a>1. Obtaining Genotype data**

You can download the effect sizes, genotype, phenotype, and summary statistics used in this simulation from [this link](https://drive.google.com/drive/folders/1MjLUdIxfneM3-Oh-5LX-AG9MzQ3PkrHL?usp=sharing).

**<a name="sim-analysis"></a>2. Running the analysis**

We will use Scene1 as the example for demonstration. You can repeat the same procedure for other Scenes.

```
cd result/EUR_EAS/Scene1

# PRS-CSx
cd PRS_CSx/; sbatch --array=1-22 PRS_CSx.sh
# after all jobs finish
sbatach PRS_CS_resx.sh

# LDpred2
cd ldpred2/; sbatch --array=1-10 ldpred.sh

# SDPRX
cd SDPRX/; sbatch --array=1-10 SDPRX.sh

# XPASS
cd XPASS/; sbatch --array=1-10 XPASS.sh
```

# <a name="real"></a>Real data applications

## <a name="real-req"></a>Requirements

* LDSC 
* PopCorn
* about 230 Gb desk space

## Workflow

**<a name="real-geno"></a>1. Obtaining Genotype data**

The list of SNPs passing quality control can be found in the directory `genotype/EAS` and `genotype/AFR`. One can then obtain the genotype data if having access to UK biobank.

**<a name="real-ref"></a>2. Constructing the reference LD matrix**

The refernce LD of PRS-CSx and SDPRX can be obtained from their website. The `ref/` directory has the script used to construct LD matrix for LDpred2.

**<a name="real-ss"></a>3. Obtaining and cleaning the summary statistics**

The study of the summary statistics is listed in the Web resources section of the manuscript. If you need assistance, feel free to submit to the issue. Here we provide an example on downloading and processing the height summary statistics from the GIANT consortium. The procedure is similar for other traits. 

```
cd UKB_real/HGT/summ_stats/

# download
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

# make sure you check the name of summary statistics is right
# you also need to change the path of LDSC, munge_summstats.py, and reference of LDSC
Rscript clean1.R
```

**<a name="real-analysis"></a>4. Running the analysis**

Here is the example on running analysis for height. The procedure is similar for other traits. 

```
cd UKB_real/HGT/result/

# XPASS
cd XPASS/; sbatch XPASS.sh

# PRS-CSX
cd ../PRS_CS/; sbatch --array=1-22 PRS_CSx.sh
# after all jobs finish
sbatch PRS_CSx_res.sh

# SDPRX
cd ../SDPR/; sbatch --array=1-22 SDPRX.sh
# after all jobs finish

# LDpred2
cd ../ldpred2/; sbatch ldpred2.sh
```

 

