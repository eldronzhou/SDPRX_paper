
library(XPASS)
library(data.table)
library(RhpcBLASctl)
library(ieugwasr)

args = commandArgs(trailingOnly=T)

j = args[1]
scene = args[2]

ref_EAS = "../../../../genotype/EAS/validate/validate_5k"

ref_EUR = "../../../../genotype/EUR/validate/validate_5k"

ss1 = paste0("../../../../summary_stat/EUR_EAS/",scene,"/EUR/discover/sim_",j,"_xpass.txt")

ss2 = paste0("../../../../summary_stat/EUR_EAS/",scene,"/EAS/discover/sim_",j,"_xpass.txt")

z1 = fread(ss1)
pval1 = data.frame(rsid=z1$SNP,pval=2*pnorm(abs(z1$Z),lower.tail=F))
clp1 = ld_clump(pval1, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
		    bfile=ref_EUR, plink_bin="plink")
snps1 = clp1$rsid

z2 = fread(ss2)
pval2 = data.frame(rsid=z2$SNP,pval=2*pnorm(abs(z2$Z),lower.tail=F))
clp2 = ld_clump(pval2, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
		                    bfile=ref_EAS, plink_bin="plink")
snps2 = clp2$rsid

fit = XPASS(file_z1=ss1, file_z2=ss2, file_ref1=ref_EUR,
		     file_ref2 = ref_EAS, snps_fe1=snps1, snps_fe2=snps2)

write.table(fit$mu[,c("SNP","A1","mu_XPASS1","mu_XPASS2")], file=paste0("sim_",j,".txt"), row.names=F, col.names=T, quote=F, append=F)

