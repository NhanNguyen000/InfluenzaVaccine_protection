# prepare data for qtl analysis based on Martijn code
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(RNOmni)
#setwd('/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters')

load('/vol/projects/CIIM/Influenza/iMED/metabolic/analysis/output/cohort_2/data.RData')
abs <- read.csv('/vol/projects/CIIM/Influenza/iMED/metadata/FCtiters_adjusted_baseline.csv', sep = ',', row.names = 1)

abs_T1 <- abs %>% select(pID, H1N1_T1, H3N2_T1, B_T1) %>% tibble::column_to_rownames('pID') %>% t()
phenotype_T1 <- apply(abs_T1, 1, function(x) {RankNorm(u = x)}) %>% t()

abs_T4 <- abs %>% select(pID, H1N1_T4, H3N2_T4, B_T4) %>% tibble::column_to_rownames('pID') %>% t()
phenotype_T4 <- apply(abs_T4, 1, function(x) {RankNorm(u = x)}) %>% t()

# Covariates
covariates <- meta %>%
  filter(time == 'T1') %>%
  select(ProbandID, age, gender) %>%
  mutate(gender = ifelse(gender == 'male', 1, 2)) %>%
  set_colnames(c('sampleID', 'age', 'gender')) %>%
  tibble::column_to_rownames('sampleID') %>%
  t()

# Genotype
dos_temp <- data.table::fread('/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/dosage.txt') %>% as.data.frame
colnames(dos_temp) <- reshape2::colsplit(colnames(dos_temp), ']', names = c('trash', 'id')) %>% pull(id)
colnames(dos_temp) <- reshape2::colsplit(colnames(dos_temp), ':', names = c('nr', 'trash')) %>% pull(nr)

dos_temp$CHROM <- substr(dos_temp$CHROM, 4, nchar(dos_temp$CHROM))
dos_temp$ID <- colsplit(dos_temp$ID, '%', names = c('hg38', 'rsID')) %>% pull(rsID)


#Hg19 rownames
rownames(dos_temp) <- paste0(dos_temp$CHROM, ':', dos_temp$POS, ':', dos_temp$ALT, ':', dos_temp$REF, '%', dos_temp$ID)
dos <- dos_temp %>% select(-CHROM, -ID, -POS, -REF, -ALT)

colnames(dos) <- lapply(strsplit(colnames(dos), split='-', fixed=T), function(x) {x[[2]]}) %>% unlist() %>% stringr::str_remove(., "^0+")

# What are the samples that we have?
samples <- intersect(colnames(phenotype_T1), colnames(dos)) #159

# Output
phenotype_T1 <- phenotype_T1[, samples] %>% as.data.frame()
phenotype_T4 <- phenotype_T4[, samples] %>% as.data.frame()
covariates <- covariates[, samples] %>% as.data.frame()
dos <- dos[, samples] %>% as.data.frame()

phenotype_T1[1:2, 1:4]
phenotype_T4[1:2, 1:4]
covariates[1:2, 1:4]
dos[1:4, 1:4]

stopifnot(identical(names(phenotype_T1), names(phenotype_T4)))
stopifnot(all(colnames(phenotype_T1) == colnames(dos)))
stopifnot(all(colnames(phenotype_T1) == colnames(dos)))
stopifnot(all(colnames(covariates) == colnames(dos)))

# Output
data.table::fwrite(file = 'processedDat/qtl/phenotypes_T1.csv', x = phenotype_T1, sep = '\t', quote = F, eol = '\n', row.names = T)
data.table::fwrite(file = 'processedDat/qtl/phenotypes_T4.csv', x = phenotype_T4, sep = '\t', quote = F, eol = '\n', row.names = T)
data.table::fwrite(file = 'processedDat/qtl/genotype.csv', x = dos, sep = '\t', quote = F, eol = '\n', row.names = T)
data.table::fwrite(file = 'processedDat/qtl/covariates.csv', x = covariates, row.names=T, quote=F, sep = '\t', eol = '\n')