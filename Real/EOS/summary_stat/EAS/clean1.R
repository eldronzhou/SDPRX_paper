library(readr)

a = read_tsv("Eosino_GWAS_in_BBJ_autosome.tsv")
d = a[,c(5,3,4,8,7)]
d$N = 62076
colnames(d) = c("SNP","A1","A2","b","p","N")

write.table(d, file="clean1.txt", row.names=F, col.names=T, quote=F, append=F)

system("module unload R; source activate ldsc; python ~/software/ldsc/munge_sumstats.py --sumstats clean1.txt --out clean; ~/software/ldsc/ldsc.py --h2 clean.sumstats.gz --ref-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eas_ldscores/ --out ./clean_h2 --w-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eas_ldscores/")

# prepare for PRS-CS
b = read.table("clean.sumstats.gz", header=T, stringsAsFactors=F)
snp = read.table("../../../snplist/1kg_prscs_inter_eas.txt", stringsAsFactors=F)
b = b[b$SNP %in% snp$V1,]
b$P = 2*pnorm(-abs(b$Z))
colnames(b) = c("SNP","A1","A2","BETA","N","P")
write.table(b[,c(1,2,3,4,6)], file="PRS_cs.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for SDPR
write.table(b[,c(1,2,3,4,6,5)], file="SDPR.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for xpass
colnames(b)[4] = "Z"
write.table(b[,c("SNP","N","Z","A1","A2")], file="XPASS.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for ldpred2
a = a[a$rsids %in% b$SNP,]
a$N = 62076
write.table(a[,c(1,5,2,4,3,8,9,7,11)], file="ldpred2.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

system("rm -rf clean1.txt")

system("module unload R; source activate python2; popcorn fit -v 0 --cfile ../../../software/Popcorn/ref/EUR_EAS_all_gen_eff.cscore --gen_effect --sfile1 ../EUR/XPASS.txt --sfile2 XPASS.txt ./popcorn.txt")

print(median(b$N))



