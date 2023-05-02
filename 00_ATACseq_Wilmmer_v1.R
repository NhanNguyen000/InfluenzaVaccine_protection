# ATACseq data
/vol/projects/CIIM/Influenza/supplemental_data/Wimmers_et_al/tiv_scatac.h5ad.gz
#ATACseq from Wimmers et al. - CHECK THE CODE: /vol/projects/BIIM/MMR_vaccine/scATACseq/ArchR/Code,
rm(list = ls())

load("/vol/projects/CIIM/Influenza/supplemental_data/Wimmers_et_al/analysis_rna.R") #  bad restore file magic number (file may be corrupted) (?)
load("/vol/projects/CIIM/Influenza/supplemental_data/Wimmers_et_al/sc_sequencing/load_data.R") #  bad restore file magic number (file may be corrupted) (?)

# read the h5ad.gz file
#BiocManager::install("rhdf5")
library("rhdf5")
h5ls("data/tiv_scatac.h5ad")
mydata <- h5read("data/tiv_scatac.h5ad", "/")
mydata <- h5read("data/tiv_scatac.h5ad", "/obs")

library(Seurat)
library(SeuratDisk)
library("anndata")
library(Seurat)
library(SeuratDisk)
library(SingleR)
library(tidyverse)
set.seed(123)
ad <- read_h5ad("data/tiv_scatac.h5ad")
dat <- Read10X_h5("data/tiv_scatac.h5ad")
f <- hdf5r::H5File$new("data/tiv_scatac.h5ad")

Convert("data/tiv_scatac.h5ad", dest = "h5seurat", overwrite = TRUE)
dat <- LoadH5Seurat("data/tiv_scatac.h5seurat") # Error: Missing required datasets 'levels' and 'values'
dat <- LoadH5Seurat("data/tiv_scatac.h5seurat", meta.data = FALSE, mics = FALSE) # also not work, Error: Missing required datasets 'levels' and 'values'
dat <- LoadH5Seurat("data/tiv_scatac.h5seurat", meta.data = FALSE, mics = FALSE, tools = FALSE) # also not work,
dat <- LoadH5Seurat("data/tiv_scatac.h5seurat", assays = "RNA")
dat <- as.Seurat("data/tiv_scatac.h5seurat")


Convert("data/tiv_scatac.h5ad", "test.h5seurat")
train_seurat <- LoadH5Seurat("test.h5seurat", assays = "RNA") # Error: Missing required datasets 'levels' and 'values'
a <- readH5AD("data/tiv_scatac.h5ad")
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
#>  2 dimensional reductions calculated: pca, umap


#Reading the recently created file into a seurat object

seu <- LoadH5Seurat("../combined/all_weeks_combined_named.h5seurat", assay = "RNA")

# read files
a <- read.table()

metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

