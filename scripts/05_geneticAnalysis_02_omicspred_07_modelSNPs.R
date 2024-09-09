# Zoom in the SNPs in OmicPred models

rm(list = ls())

library(tidyverse)

# load & prepare data ===========================================
qtl_dat <- read.table("/vol/projects/CIIM/Influenza/ZirrFlu/genetics_abqtl/out/mapping/H3N2_T4.txt", 
                header = TRUE)

#  CCL7 in somalogic - OPGS000414 model
folder_path <- "/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/omicsPred_models/"
CCL7_model <- read.table(paste0(folder_path, "SomaScan/OPGS000414_model.txt"), header = TRUE)


CCL7_snps <- c()
for (i in 1:length(CCL7_model$rsid)) {
  CCL7_snps <- c(CCL7_snps, grep(CCL7_model$rsid[1], qtl_dat$SNP))
}

CCL7snps_qtlDat <- qtl_dat[CCL7_snps, ] # No pvalue < 0.05

