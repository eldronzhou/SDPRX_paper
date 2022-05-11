
library(stringr)
library(data.table)
library(bigsnpr)
library(methods)

bigparallelr::set_blas_ncores(1)
options(bigstatsr.check.parallel.blas = FALSE)

ss = fread("../../summary_stat/EAS/ldpred2.txt")

obj.bigsnp = snp_attach("../../../ref/EAS/1000G_phase3_common_norel.rds")
G = obj.bigsnp$genotypes
CHR = as.integer(obj.bigsnp$map$chromosome)

colnames(ss) = c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
map = obj.bigsnp$map[-c(2:3)]
names(map) = c("chr", "pos", "a0", "a1")
info_snp = snp_match(ss, map, strand_flip=FALSE)

tmp = tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

for (chr in 1:22) {
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]
    ind.chr3 = match(ind.chr2, which(CHR == chr))

    # read corr
    corr0 = readRDS(paste0("../../../ref/EAS/ldpred2/1kg_mh3_ldpred2_chr",chr,".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
	df_beta = info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
	ld = Matrix::colSums(corr0^2)
	corr = as_SFBM(corr0, tmp)
    }
    else {
	df_beta = rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
	ld = c(ld, Matrix::colSums(corr0^2))
	corr$add_columns(corr0, nrow(corr))
    }
} 

(ldsc = with(df_beta, snp_ldsc(ld, length(ld), chi2=(beta/beta_se)^2, 
			       sample_size = n_eff, blocks = NULL)))

h2_est <- ldsc[["h2"]]

# LDpred2-inf
beta_inf = snp_ldpred2_inf(corr, df_beta, h2_est)
res = data.frame(SNP=info_snp$rsid, A1=info_snp$a1,
                 BETA=beta_inf)

# LDpred2-grid
(h2_seq = round(h2_est * c(0.7, 1, 1.4), 4))
(p_seq = signif(seq_log(1e-5, 1, length.out=21), 2))
(params = expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE)))

beta_grid = snp_ldpred2_grid(corr, df_beta, params, ncores=10)

#for (p in 1:ncol(beta_grid)) {
#    res = data.frame(SNP=info_snp$rsid, A1=info_snp$a1, BETA=beta_grid[,p])
#    write.table(res, file=paste0("ldpred_",p,".txt"), row.names=F, col.names=F, quote=F, append=F)
#}

# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init=h2_est,
	vec_p_init = seq_log(1e-4, 0.9, 10), ncores=10)

beta_auto = sapply(multi_auto, function(auto) auto$beta_est)
pred_auto = big_prodMat(G, beta_auto, ind.col=df_beta[["_NUM_ID_"]])
sc = apply(pred_auto, 2, sd)
keep = abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto = rowMeans(beta_auto[, keep])

res = data.frame(SNP=info_snp$rsid, A1=info_snp$a1, BETA=cbind(beta_inf, beta_grid, final_beta_auto))
res[is.na(res)] = 0
write.table(res, file=paste0("ldpred.txt"), row.names=F, col.names=F, quote=F, append=F)



