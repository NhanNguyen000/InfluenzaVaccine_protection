# run qtl analysis based on Martijn code, 
# but after Maritjn talk with Javi for the genetic imputation, he re-run the code with Javi pipeline
# no the correct QTL result is in /vol/projects/CIIM/Influenza/ZirrFlu/genetics_abqtl/out/mapping/, and only for T4

rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(MatrixEQTL)

# load data =================================================
SNP_fl = 'processedDat/qtl/genotype.csv'
covariates_fl = 'processedDat/qtl/covariates.csv'

# setup parameter ----------------------------------
p.val = 1

useModel <- modelLINEAR
errorCovariance <- numeric()

# snp settings
snps = SlicedData$new()
snps$fileOmitCharacters = "NA"
snps$fileDelimiter = "\t"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 200

# covariates settings
cvrt = SlicedData$new()
cvrt$fileOmitCharacters = "NA"
cvrt$fileDelimiter = "\t"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 200

# phenotyple settings
pheno = SlicedData$new()
pheno$fileOmitCharacters = "NA"
pheno$fileDelimiter = "\t"
pheno$fileSkipRows = 1
pheno$fileSkipColumns = 1
pheno$fileSliceSize = 200

# run qtl analysis -------------------------
for (time in c("T1", "T4")) {
  pheno_fl = paste0("processedDat/qtl/phenotypes_", time, ".csv")
  output_fl = paste0("processedDat/qtl/outcome_", time, ".csv")
  
  # add input to the settings
  snps$LoadFile(SNP_fl)
  cvrt$LoadFile(covariates_fl)
  pheno$LoadFile(pheno_fl)
  
  # run qtl
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
  
}




