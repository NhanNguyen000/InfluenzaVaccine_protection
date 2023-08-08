rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
# scRNAseq data (Dusseldorf cohort) ------------------------------------------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData") # Note: marker.list is based on accumulated previous information (no specific papers)
metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

unique(data$vaccination)
unique(data$celltype)

before_vac <- subset(x = data, subset = vaccination == "before") # subset only scRNAseq before vaccination
after_vac <- subset(x = data, subset = vaccination == "after") # subset only scRNAseq before vaccination

# specific protein ---------------------
protein <- "CD83"
before_vac %>% 
  VlnPlot(features =protein, group.by = "celltype") + NoLegend()

after_vac %>% 
  VlnPlot(features =protein, group.by = "celltype") + NoLegend()

data %>%
  VlnPlot(features =protein, split.by = "vaccination", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)

# UMAP --------------------------------------
before_vac %>% 
  DimPlot(group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()

before_vac %>% 
  FeaturePlot(reduction = "umap", features = protein, # specific protein
            sort.cell = TRUE, min.cutoff = '1', label = TRUE) 



