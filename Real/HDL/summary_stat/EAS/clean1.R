a = readr::read_tsv("HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")
snp = read.table("../../../snplist/1kg_prscs_inter_eas.txt")
a[a$pvalue<1e-323,"pvalue"] = 1e-323
a = a[a$rsID %in% snp$V1, ]
d = a[, c("rsID","REF","ALT","EFFECT_SIZE","pvalue","N")]
d = d[!duplicated(d$rsID),]
colnames(d) = c("SNP","A1","A2","beta","P","N")
write.table(d, file="clean1.txt", row.names=F, col.names=T, quote=F, append=F)

system("module unload R; source activate ldsc; python ~/software/ldsc/munge_sumstats.py --sumstats clean1.txt --out clean")

system("module unload R; source activate ldsc; ~/software/ldsc/ldsc.py --h2 clean.sumstats.gz --ref-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/ --out ./clean_h2 --w-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/")

# prepare for PRS-CS
b = read.table("clean.sumstats.gz", header=T, stringsAsFactors=F)
d = d[d$SNP %in% b$SNP,]
stopifnot(all.equal(d$SNP, b$SNP))
d$beta = b$Z
colnames(d) = c("SNP","A1","A2","BETA","P","N")
write.table(d[,c(1,2,3,4,5)], file="PRS_cs.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for SDPR
write.table(d, file="SDPR.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for XPASS
colnames(d)[4] = "Z"
write.table(d[,c("SNP","N","Z","A1","A2")], file="XPASS.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for ldpred2
d = a[a$rsID %in% b[,1], c("CHROM","rsID","POS_b37","REF","ALT","EFFECT_SIZE","SE","pvalue","N")]
write.table(d, file="ldpred2.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)
system("rm -rf clean1.txt")

system("module unload R; source activate python2; popcorn fit -v 0 --cfile ../../../software/Popcorn/ref/EUR_EAS_all_gen_eff.cscore --gen_effect --sfile1 ../EUR/XPASS.txt --sfile2 XPASS.txt ./popcorn.txt")

print(median(d$N))



