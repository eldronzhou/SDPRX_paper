
library(stringr)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
k = args[3] # scene

pops = c("EUR","EAS")

for (pop in pops) {

pheno = read.table(paste0("../../../../phenotype/EUR_EAS/",k,"/",pop,"/validate/sim_",j,".phen"), header=F)

file = list.files(path=paste0("sim_",j,"/",pop,"/"), pattern="*.txt$", full.names=T)
r2 = rep(0, length(file))

for (i in seq_along(file)) {
	command = paste0("module load PLINK/1.90-beta5.3; plink --bfile ../../../../genotype/",pop,"/validate/validate_5k --score ",file[i]," 1 2 3 --out ",file[i] )
	system(command, intern=F, wait=T)
	dat = read.table(paste0(file[i],".profile"), header=T)
	r2[i] = cor(pheno[,3], dat$SCORE)^2
}

idx = which.max(r2)
print(paste0(file[idx], ": ", r2[idx]))


#name = str_extract(file[idx], "r2_0.[0-9]*")

command = paste0("module load PLINK/1.90-beta5.3; plink --bfile ../../../../genotype/",pop,"/test/test_5k --score ", file[idx], " 1 2 3  --out ../predict/ldpred2/sim_",j,"_",pop)

system(command, wait=T, intern=F)
}

