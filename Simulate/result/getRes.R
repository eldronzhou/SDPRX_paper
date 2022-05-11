library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)


path <- "/home/tc732/project/SDPRX/"
path_sim <- paste0(path, "simulation_1000GP_Phase3_chr1/")

method <- c("SDPR", "PRS_CSx")
res <- data.frame(pops = character(0), 
                  method = character(0),
                  pop = character(0),
                  R2 = numeric(0))

pops_ls <- list(c("EUR", "EAS"), c("EUR", "AFR"))

for(j in 1:length(pops_ls)){
  pop1 <- pops_ls[[j]][1]
  pop2 <- pops_ls[[j]][2]
  pops <- paste0(pop1, "_", pop2)
  
  for (i in 1:10) {
    pheno_pop1 <- read.table(paste0(path_sim, "phenotype/", pops, "/Scene1/", 
                                    pop1, "/test/sim_",i,".phen"), 
                             header=F, stringsAsFactors=F)
    pheno_pop2 <- read.table(paste0(path_sim, "phenotype/", pops, "/Scene1/", 
                                    pop2, "/test/sim_",i,".phen"), 
                             header=F, stringsAsFactors=F)
    # SDPR  <- read.table(paste0(sz[j],"/predict/SDPR/sim_",i,".profile"), header=T, stringsAsFactors=F)
    # stopifnot(all.equal(SDPR$FID, pheno$V1))
    # res[50*(j-1)+5*(i-1)+1,1] = sz[j]
    # res[50*(j-1)+5*(i-1)+1,2] = method[1]
    # res[50*(j-1)+5*(i-1)+1,3] = cor(pheno$V3, SDPR$SCORE)^2
    
    prs_csx_pop1 <- read.table(paste0(path_sim,"/result/", pops ,"/Scene1/PRS_CSx/", 
                                      pop1 ,"/sim_",i,".profile"), 
                               header=T, stringsAsFactors=F)
    prs_csx_pop2 <- read.table(paste0(path_sim,"/result/", pops ,"/Scene1/PRS_CSx/", 
                                      pop2 ,"/sim_",i,".profile"), 
                               header=T, stringsAsFactors=F)
    
    stopifnot(all.equal(prs_csx_pop1$FID, pheno_pop1$V1))
    res <- rbind(res, data.frame(pops=pops, 
                                 method="PRS_CSx", 
                                 pop=pop1, 
                                 R2=cor(pheno_pop1$V3, prs_csx_pop1$SCORE)^2))
    res <- rbind(res, data.frame(pops=pops, 
                                 method="PRS_CSx", 
                                 pop=pop2, 
                                 R2=cor(pheno_pop2$V3, prs_csx_pop2$SCORE)^2))
    
  }
  
}


# tiff("Scenario_1A.tiff", units="in", width=6, height=4, res=300)
res <- res %>% group_by(pops, method, pop) %>% mutate(med=median(R2))
ggplot(data=res, aes(x=pop, y=R2,fill=method)) + 
  geom_boxplot()
  # geom_boxplot(outlier.shape=NA, position=position_dodge(width=.75), width=.4) +
  # geom_point(position=position_jitterdodge(dodge.width=.75), color="black", size=.1) +
  # ggtitle("Scenario 1A") + theme(plot.title=element_text(hjust = 0.5),
  #                                text=element_text(size=10))
# dev.off()
