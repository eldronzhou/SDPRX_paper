
library(methods)
library(bigsnpr)
library(stringr)

args = commandArgs(trailingOnly = T)

i = args[1] # chr

obj.bigSNP = snp_attach("../1000G_phase3_common_norel.rds")

G = obj.bigSNP$genotypes
CHR = obj.bigSNP$map$chromosome
POS = obj.bigSNP$map$physical.pos

POS2 = snp_asGeneticPos(CHR, POS)

idx = which(CHR == i)
corr = snp_cor(G, ind.col=idx, infos.pos=POS2[idx], size=3/1000)

saveRDS(corr, file=paste0("1kg_mh3_ldpred2_chr",i,".rds"))

