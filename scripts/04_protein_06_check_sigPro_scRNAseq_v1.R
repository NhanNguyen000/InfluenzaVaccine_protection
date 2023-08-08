rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
# scRNAseq data (Dusseldorf cohort) ------------------------------------------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData") # Note: marker.list is based on accumulated previous information (no specific papers)
metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

data <- RegroupIdents(data, metadata = "celltype") # rename cluster by cell type
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


# dotplot for all related protein to CD83 ---------------------------
CD83_relatedPro <- c("CD83", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD40", "CD80", "CD86", "TNF", "CCR7") # based on stringDB

DotPlot(object = before_vac, features = CD83_relatedPro) + 
  theme(axis.text.x = element_text(angle = 90))
