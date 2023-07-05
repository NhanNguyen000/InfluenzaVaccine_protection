try(dev.off())
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
abs %<>% select(pID, H1N1_T4, H3N2_T4, B_T4) %>% tibble::column_to_rownames('pID') %>% t()
phenotype <- apply(abs, 1, function(x) {RankNorm(u = x)})
phenotype <- t(phenotype)

# Covariates
covariates <- meta %>%
  filter(time == 'T1') %>%
  select(ProbandID, age, gender) %>%
  mutate(gender = ifelse(gender == 'male', 1, 2)) %>%
  set_colnames(c('sampleID', 'age', 'gender')) %>%
  tibble::column_to_rownames('sampleID') %>%
  t()

# Genotype
data.table::fread('/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/dosage.txt') %>% as.data.frame -> dos
colnames(dos) <- reshape2::colsplit(colnames(dos), ']', names = c('trash', 'id')) %>% pull(id)
colnames(dos) <- reshape2::colsplit(colnames(dos), ':', names = c('nr', 'trash')) %>% pull(nr)

dos$CHROM <- substr(dos$CHROM, 4, nchar(dos$CHROM))
dos$ID <- colsplit(dos$ID, '%', names = c('hg38', 'rsID')) %>% pull(rsID)


#Hg19 rownames
rownames(dos) <- paste0(dos$CHROM, ':', dos$POS, ':', dos$ALT, ':', dos$REF, '%', dos$ID)
dos %<>% select(-CHROM, -ID, -POS, -REF, -ALT)

colnames(dos) <- lapply(strsplit(colnames(dos), split='-', fixed=T), function(x) {x[[2]]}) %>% unlist() %>% stringr::str_remove(., "^0+")

# What are the samples that we have?
intersect(colnames(phenotype), colnames(dos)) -> samples #159

# Output
phenotype <- phenotype[, samples] %>% as.data.frame()
covariates <- covariates[, samples] %>% as.data.frame()
dos <- dos[, samples] %>% as.data.frame()

covariates[1:2, 1:4]
phenotype[1:2, 1:4]
dos[1:4, 1:4]

stopifnot(all(colnames(phenotype) == colnames(dos)))
stopifnot(all(colnames(covariates) == colnames(dos)))

# Output
data.table::fwrite(file = 'phenotypes.csv', x = phenotype, sep = '\t', quote = F, eol = '\n', row.names = T)
data.table::fwrite(file = 'genotype.csv', x = dos, sep = '\t', quote = F, eol = '\n', row.names = T)
data.table::fwrite(file = 'covariates.csv', x = covariates, row.names=T, quote=F, sep = '\t', eol = '\n')