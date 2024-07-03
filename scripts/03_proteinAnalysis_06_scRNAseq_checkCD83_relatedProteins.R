rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
# scRNAseq data (Dusseldorf cohort) ------------------------------------------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData") # Note: marker.list is based on accumulated previous information (no specific papers)
metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

data <- RegroupIdents(data, metadata = "celltype") # rename cluster by cell type
unique(data$vaccination) # check time points
unique(data$celltype) # check celly types

before_vac <- subset(x = data, subset = vaccination == "before") # subset only scRNAseq before vaccination
after_vac <- subset(x = data, subset = vaccination == "after") # subset only scRNAseq before vaccination

# overview about cell types
before_vac %>% 
  DimPlot(group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()

# plot for specitic protein ----------------------------------------------------
protein <- "CD83"
protein <- "CD40"
protein <- "CD80"
protein <- "CCR7"

# UMAP
umap_protein <- before_vac %>% 
  FeaturePlot(reduction = "umap", features = protein, # specific protein
              order = TRUE, min.cutoff = '1', label = TRUE, label.size = 5)
umap_protein

# violine plot per cell type
before_vac %>% 
  VlnPlot(features =protein, group.by = "celltype") + NoLegend()

after_vac %>% 
  VlnPlot(features =protein, group.by = "celltype") + NoLegend()

data %>%
  VlnPlot(features =protein, split.by = "vaccination", group.by = "celltype", 
          combine = FALSE)


# dotplot for all related protein to CD83 ---------------------------
CD83_relatedPro <- c("CD83", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD40", "CD80", "CD86", "TNF", "CCR7") # based on stringDB

dotplot_CD83_relatedPro <- DotPlot(object = before_vac, features = CD83_relatedPro) + 
  theme( axis.text.x = element_text(angle = 90, size = 16),
         axis.text.y = element_text(size = 16))

dotplot_CD83_relatedPro

# save the plot  -------------------------------------------------
ggsave(paste0("output/umap_scRNAseq_", protein, ".png"), 
       umap_protein, device = "png")

ggsave("output/dotplot_scRNAseq_CD83_relatedProteins.png", 
       dotplot_CD83_relatedPro, device = "png")
