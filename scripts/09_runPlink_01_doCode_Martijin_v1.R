rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(MatrixEQTL)

setwd('/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters')

p.val = 1

SNP_fl = 'genotype.csv'
covariates_fl = 'covariates.csv'
pheno_fl = 'phenotypes.csv'
output_fl = 'out.csv'

useModel <- modelLINEAR
errorCovariance <- numeric()

snps = SlicedData$new()
snps$fileOmitCharacters = "NA"
snps$fileDelimiter = "\t"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 200
snps$LoadFile(SNP_fl)

cvrt = SlicedData$new()
cvrt$fileOmitCharacters = "NA"
cvrt$fileDelimiter = "\t"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 200
cvrt$LoadFile(covariates_fl)

pheno = SlicedData$new()
pheno$fileOmitCharacters = "NA"
pheno$fileDelimiter = "\t"
pheno$fileSkipRows = 1
pheno$fileSkipColumns = 1
pheno$fileSliceSize = 200
pheno$LoadFile(pheno_fl)

quest_qtl = Matrix_eQTL_engine(
  snps = snps,
  gene = pheno,
  cvrt = cvrt,
  output_file_name = output_fl,
  pvOutputThreshold = p.val,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = T,
  min.pv.by.genesnp = F,
  noFDRsaveMemory = T,
  pvalue.hist = T
)
