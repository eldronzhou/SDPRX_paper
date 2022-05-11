
library(XPASS)
library(data.table)
library(RhpcBLASctl)
library(ieugwasr)


ref_EAS = "../../../ref/EAS/1000G_phase3_common_norel"

ref_EUR = "/ysm-gpfs/pi/zhao/gz222/1000g_phase3/genotype_1KG_eur_SNPmaf5/hm3/eur_SNPmaf5_nomhc"

ss1 = paste0("../../summary_stat/EAS/XPASS.txt")

ss2 = paste0("../../summary_stat/EUR/XPASS.txt")

z1 = fread(ss1)
pval1 = data.frame(rsid=z1$SNP,pval=2*pnorm(abs(z1$Z),lower.tail=F))
clp1 = ld_clump(pval1, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
		    bfile=ref_EAS, plink_bin="plink")
snps1 = clp1$rsid

z2 = fread(ss2)
pval2 = data.frame(rsid=z2$SNP,pval=2*pnorm(abs(z2$Z),lower.tail=F))
clp2 = ld_clump(pval2, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
		                    bfile=ref_EUR, plink_bin="plink")
snps2 = clp2$rsid

fit = XPASS(file_z1=ss1, file_z2=ss2, file_ref1=ref_EAS,
		     file_ref2=ref_EUR, snps_fe1=snps1, snps_fe2=snps2)

write.table(fit$mu[,c("SNP","A1","mu_XPASS1","mu_XPASS2")], file="xpass.txt", row.names=F, col.names=T, quote=F, append=F)

