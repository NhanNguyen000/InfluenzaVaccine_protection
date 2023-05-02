# Readme ----------------------------------------------------------
Created on: Thu May 20 15:47:25 CEST 2021
Last updated on: May 20th

Created by: Martijn Zoodsma (Helmholtz Centre for Infection Research, Hannover, Germany)
Principal investigator: Prof. Dr. Yang Li (Helmholtz Centre for Infection Research. Hannover, Germany)
Contact information: martijn.zoodsma@helmholtz-hzi.de

Metadata: ./data/metadata.csv
Raw data: /vol/data/gmak/runs/raw/novaseq/210330_A00278_0261_BH2YY7DRXY

# Folders
/vol/projects/mzoodsma/influenza/scRNAseq/analysis
└── output
├── cell_proportions
├── celltypes
├── DEG
└── enrichment

# Software
R

# Description
Contains R scripts used for scRNAseq analysis. Scripts run.sh and run_all.sh may be used to quickly submit R scripts to the SGE cluster.

Output is stored in ./output/. Contains:
- UMAPs from Seurat
- Differential expression results and visualization (volcano plots)
- Enrichment analysis full results and visualization
- Cell proportion analysis
- Raw Seurat object is stored in : data.RData for further use.

R scripts:
- 1. load_data.R: Load the raw data, produce Seurat object and UMAPs
- 2. Annotation: Celltype annotation of Seurat clusters
- 3. DE.R: Differential expression analysis
- 4. cell_proportions.R: Dirichlet regression to check for differences in cell proportion over conditions
- 5. enrichment.R: GO and KEGG enrichment for the differentially expressed genes within a celltype.

# 1. load_data.R ----------------------------------------------------------
#Loading data and merging into Seurat object 

rm(list = ls())
try(dev.off())

set.seed(42)

options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

# Loading libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(future))

plan(strategy = "multicore", workers = 12)

OUTPUT_DIR = "/vol/projects/mzoodsma/influenza/scRNAseq/analysis/output/"
DATA_DIR = "/vol/projects/mzoodsma/influenza/scRNAseq/data/"

# Replace later on with souporcell genotype clusters
sampleIDs <- list(
  c("0" = "CP0027b", "1" = "CP0043b", "2" = "CP007",   "3" = "CP009",   "4" = "CP0028b"),
  c("0" = "CP009b",  "1" = "CP0011",  "2" = "CP0014",  "3" = "CP007b",  "4" = "CP0013"),
  c("0" = "CP0013b", "1" = "CP0018",  "2" = "CP0022",  "3" = "CP0014b", "4" = "CP0011b"),
  c("0" = "CP0028",  "1" = "CP0022b", "2" = "CP0018b", "3" = "CP0027",  "4" = "CP0043")
)

seur_list <- list()

for (i in 1:length(sampleIDs)) {
  print(paste0('Pool: ', i))
  
  # Reading in cluster assignments
  clusters <- read.table(paste0(DATA_DIR, 'pool', i, '_clusters.tsv'), header=T, stringsAsFactors = F)
  clusters$barcode <- substr(clusters$barcode, start = 1, stop = nchar(clusters$barcode)-2)
  print(table(clusters$status))
  # Reading in 10X data
  seur <- Read10X(paste0(DATA_DIR, 'pool', i, '_filtered_feature_bc_matrix'))
  colnames(seur) <- substring(colnames(seur),1,nchar(colnames(seur))-2)
  
  # Subset to take only singletons
  singletons <- subset(clusters, status == 'singlet')
  seur <- seur[, colnames(seur) %in% singletons$barcode]
  
  # Calculate % of mitochondrial genes, and remove them
  MT.index <- grep(pattern = "^MT-", x = rownames(seur), value = FALSE)
  all.sum <- Matrix::colSums(seur)
  percent.MT <- Matrix::colSums(seur[MT.index, ])/all.sum
  seur <- seur[-MT.index, ]
  
  # Calculate % of ribosomal genes, and remove them
  RG.index <- grep("^RP[SL][0-9]+$", x = rownames(seur), perl=TRUE, value = FALSE)
  percent.RG <- Matrix::colSums(seur[RG.index, ])/all.sum
  seur <- seur[-RG.index, ]
  
  # Adding the cluster assignment
  singletons$assignment <- as.factor(singletons$assignment)
  samples <- sampleIDs[[i]]
  singletons$assignment <- as.character(revalue(singletons$assignment, samples))
  
  # Basic sanity checks
  stopifnot(all(singletons$barcode == names(percent.MT)))
  stopifnot(all(singletons$barcode == names(percent.RG)))
  stopifnot(all(names(percent.MT) == names(percent.RG)))
  
  seur <- CreateSeuratObject(seur,
                             min.cells = 5,
                             meta.data = data.frame(percent.mt = percent.MT,
                                                    percent.rg = percent.RG,
                                                    sample = singletons$assignment,
                                                    status = singletons$status),
                             project = "10X")
  
  seur$pool <- i
  seur_list[[i]] <- seur
}

# Merge the seurat objects
data <- merge(seur_list[[1]],
              y = c(seur_list[[2]], seur_list[[3]], seur_list[[4]]),
              add.cell.ids = c("pool_1", "pool_2", "pool_3", "pool_4"),
              project = '10X'
)

##### REMOVE POOL 3 immediately
print(dim(data))
data <- subset(data, pool != 3)
print(table(data@meta.data$pool))
print(dim(data))
# TODO: Add the rest of the metadata to the Seurat object

data@meta.data$vaccination <- NA
meta <- read.csv(paste0(DATA_DIR, 'metadata.csv'), sep = ';', header=T)

for (sample in meta$sampleID) {
  data@meta.data[which(data@meta.data$sample == sample), 'vaccination'] <- meta[which(meta$sampleID == sample), 'vaccination']
}


Idents(data) <- 'pool'

saveRDS(object = data,
        file = paste0(OUTPUT_DIR, 'merged_seurat_raw.rds'))

# QC metrics before QC filtering
pdf(paste0(OUTPUT_DIR, 'qcmetrics_raw.pdf'), width = 10, height = 5)
print(VlnPlot(data,
              pt.size = 0.05,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rg"),
              ncol = 2))
dev.off()
#
# # Subsetting based on the  QC plot made above.
data <- subset(data,
               nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 0.25 & nCount_RNA < 5000)

pdf(paste0(OUTPUT_DIR, 'qcmetrics_filtered.pdf'), width = 10, height = 5)
print(VlnPlot(data,
              pt.size = 0.05,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rg"),
              ncol = 2))
dev.off()

pdf(paste0(OUTPUT_DIR, 'qcmetrics_filtered_split_sample.pdf'), width = 10, height = 5)
print(VlnPlot(data,
              pt.size = 0.05,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rg"),
              ncol = 2,
              split.by = 'sample'))
dev.off()

data <- NormalizeData(data, normalization.method = 'LogNormalize')
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures=2000)
data <- ScaleData(data, features = VariableFeatures(data), vars.to.regres = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))
data <- RunPCA(data, dims = 1:20)
data <- RunUMAP(data, dims = 1:20)
data <- FindNeighbors(data, reduction = 'pca', dims = 1:20)
data <- FindClusters(data, resolution = 0.65)

# PCA elbow and variable feature plots
pdf(paste0(OUTPUT_DIR, 'elbowplot_pca.pdf'))
print(ElbowPlot(data, reduction = 'pca', ndims = 20))
dev.off()

pdf(paste0(OUTPUT_DIR, 'variable_features.pdf'))
print(VariableFeaturePlot(data))
dev.off()

# UMAPs for samples, pools and Seurat clusters
pdf(paste0(OUTPUT_DIR, "sample_UMAP.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "sample", pt.size = 0.5, reduction = "umap") +
        ggtitle(label = "UMAP by samples") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf(paste0(OUTPUT_DIR, "pool_UMAP.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "pool", pt.size = 0.5, reduction = "umap") +
        ggtitle(label = "UMAP by pool") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf(paste0(OUTPUT_DIR, "pool_split_UMAP.pdf"), width = 24, height = 14)
print(DimPlot(data, split.by = "pool", pt.size = 0.5, reduction = "umap", ncol=2) +
        ggtitle(label = "UMAP by pool") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf(paste0(OUTPUT_DIR, "sample_split_UMAP.pdf"), width = 20, height = 20)
print(DimPlot(data, split.by = "sample", pt.size = 0.5, reduction = "umap", ncol=2) +
        ggtitle(label = "UMAP by pool") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf(paste0(OUTPUT_DIR, "seurat_clusters_UMAP.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "seurat_clusters", pt.size = 0.5, reduction = "umap", label = T) +
        NoLegend() +
        ggtitle(label = "UMAP by seurat clusters") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()


pdf(paste0(OUTPUT_DIR, "condition_UMAP.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "vaccination", pt.size = 0.5, reduction = "umap") +
        NoLegend() +
        ggtitle(label = "UMAP by vaccination status") +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()


# Save only the Seurat object
rm(clusters, seur, seur_list, singletons, all.sum, i, MT.index, percent.MT, percent.RG, RG.index, samples)
save.image(file = paste0(OUTPUT_DIR, 'data.RData'))

# 2. annotation.R -------------------------------------
rm(list = ls())
try(dev.off())

library(Seurat)
library(ggplot2)
library(gridExtra)
library(future)
library(dplyr)

load('/vol/projects/mzoodsma/influenza/scRNAseq/analysis/output/data.RData')
plan(strategy = "multicore", workers = 12)


## Multiple steps:
# 1. Generate dotplots for the PBMC markers
# 2. Generate dotplots for the top 20 markers per celltype
# 3. Map celltypes to clusters based on plots generate in 1. and 2.

# 1. Dotplots for PBMC markers
# Marker Genes
# We use the following list of marker genes to annotate our cells
# Source: Suppressive myeloid cells are a hallmark of severe COVID-19, Cell 2020
# IL7R,CCR7,TCF7  Naive CD4 T
# IL7R            Memory CD4 T
# CD8A,GZMK       CD8 T
# NKG7,GZMB       NK
# CD14, LYZ       CD14+ monocytes
# FCGR3A          CD16+ monocytes
# CST3, CD86      MDC, DC
# IL1B            proinflammatory
# CD79A           B cell
# CD27,SDC1       plasma
# PPBP            Mk
# KIT,TPSAB1      Mast

#PBMC narkers
# Source: obmc3k vigenette seurat
# Markers 	      Cell Type
# IL7R, CCR7 	    Naive CD4+ T
# IL7R, S100A4 	  Memory CD4+
# CD14, LYZ 	    CD14+ Mono
# MS4A1 	        B
# CD8A 	          CD8+ T
# FCGR3A, MS4A7 	FCGR3A+ Mono
# GNLY, NKG7 	    NK
# FCER1A, CST3 	  DC
# PPBP 	          Platelet

#some gene markers for general PBMC annotation,
#IL7R+,CCR7+,TCF7+                            Naive CD4 T
#IL7R,CCR7-,TCF7-                             Memory CD4 T
#CD8A+,CD8B+,GZMK+                            CD8 T
#NKG7+,NCAM1+,CD8A-                           NK
#CD14+, LYZ+, FCGR3A-                         CD14+ monocytes
#FCGR3A+, CD14-                               CD16+ monocytes
#CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-       mDC, pDC (dendritic cells)
#IL1B,TNF,IL6                                 proinflammatory markers  # macrophages
#CD79A+, CD27-, SDC1-                         B cell
#CD79A+, CD27+, SDC1+                         plasmablast
#PPBP+                                        megakaryocyte/platelet
#KIT+, TPSAB1+                                Mast

# Pre-Neutrophiles, Immature Neutrophiles and Mature Neutrophiles all express the following genes
# S100A8, S100A9
# Pre-Neutrophiles express the following genes
# ELANE, MPO, PRTN3
# Pre-Neutrophiles and Immature Neutrophiles express the following genes
# PADI4, CD24, ANXA1, BPI, LCN2, ARG1, GSN, OLFM4
# Mature Neutrophiles do not express most of the marker genes. They do express S100A8 and S100A9.
# Mature Neutrophiles express the following genes, but only at low levels.
# PADI4, ANXA1, ARG1, GSN
# CD274 is only lowly expressed by mature neutrophiles.
# Source: Figure 4, Suppressive myeloid cells are a hallmark of severe COVID-19, Cell 2020


#MNC Population	Markers expressed (positive) or not expressed (negative)
#                                                  CD3
#                                                  CD19
#                                                  CD56
#                                                  CD10
#                                                  CD14
#                                                  CD66b
#                                                  CD335
#                                                  CD11c	  CD34	      CD38	   CD45RA   	CD90    	CD10	    CD123   	  CD110   	CD135
#HSC  Hematopoietic Stem Cell	                    negative  positive	negative	negative	positive	negative      -           -         -
#MPP	Multipotent Progenitor	                    negative  positive	negative	negative	negative	negative      -           -         -
#LMPP	Lymphoid-Primed Multipotent Progenitor	    negative  positive	negative	positive	negative	negative      -           -         -
#MLP	Multi Lymphoid Progenitor	                  negative  positive	negative	positive	negative	positive	    -           -         -
#CMP	Common Myeloid Progenitor	                  negative  positive	positive	negative	negative	negative	intermediate	  -	      positive
#MEP	Megakaryocyte-Erythrocyte Progenitor	      negative  positive	positive	negative	negative	negative	negative	    positive	negative
#GMP	Granulocyte-Macrophage Progenitor	          negative	positive  positive	positive	negative	negative	intermediate	  -	      positive
#CLP	Common Lymphoid Progenitor	                negative  positive	positive	positive	negative	positive	    -           -         -

# celltype_marker <- c("IL7R","CCR7","TCF7","CD8A","GZMK","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3","CD86","IL1B","CD79A","CD27","SDC1","PPBP","KIT","TPSAB1",
#                      "S100A8","S100A9","ELANE","MPO","PRTN3","PADI4","CD24","ANXA1","BPI","LCN2","ARG1","GSN","OLFM4", "CD274")
markers_pbmc <- c("CD34","CD45","IL7R","CCR7","S100A4","TCF7","CD8A","GZMK","GNLY","NKG7","GZMB","CD14","LYZ","MS4A1","MS4A7","FCER1A","FCGR3A","CST3","CD86","IL1B","CD79A","CD27","SDC1","PPBP","KIT","TPSAB1")
markers_neutrophil <- c("S100A8","S100A9","ELANE","MPO","PRTN3","PADI4","CD24","ANXA1","BPI","LCN2","ARG1","GSN","OLFM4", "CD274")
markers_bmmc <- c("CD45", "CD19", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138", "CD11c", "CD123", "CD3", "CD14", "CD16","HLA-DR", "CD141", "CD1c", "CD34", "CD117", "CD56", "CD33", "CD13", "CD11B", "CD64", "CD7", "CD8", "CD45RA", "CD69", "CD103", "CD4")
markers_mnc <- c("CD3", "CD19", "CD56", "CD10", "CD14", "CD66b", "CD335", "CD11c", "CD34", "CD38", "CD45RA", "CD90", "CD123", "CD110", "CD135")

markers_covid19 <- c("CCR7", "LEHF1", "FHIT", "PIK3IP1", "LDHB", "MAL", "ADTRP", "EEF1B2", "RPS3A", "TCF7", "LTB", "IL7R", "AQP3", "JUNB", "TRADD", "RPSA", "RPLP0", "RPS6", #CD4 T 18 markers
                     "VCAN", "LYZ", "S100A9", "S100A12", "LGALS1", "MAFB", "CD14", "PLBD1", "FCN1", "S100A8", #Classical Monocytes 10 markers (total until here 27)
                     "GZMH", "CCL5", "NKG7", "FGFBP2", "CD8A", "GZMA", "GNLY", "CST7", "CD8B", "CCL4", # CD8 T 10 markers (total until here 38)
                     "GZMB", "PRF1", "SPON2", "CLIC3", "CTSW", "HOPX", # NK cells 6 markers (total until here 44)
                     "IFII27", "IFITM3", "APOBEC3A", "TYMP", "ISG15", "CST3", "LGALS2", "IFI6", # Inflammatory Monocytes 8 markers (total until here 52)
                     "CD79A", "MS4A1", "CD74", "CD79B", "HLA-DRA", "TCL1A", "HLA-DQB1", "HLA-DPA1", "HLA-DQA1", "HLA-DPB1", # B cells 10 markers (total until here 62)
                     "GZMK", "DUSP2", "LYAR", "KLRG1", "IL32", "KLRB1", "S100B", "RPS5", "NELL2", "NOSIP", # CD 8 T 10 markers (total until here 72)
                     "FCGR3B", "NAMPT", "IFITM2", "CXCL8", "NEAT1", "G0S2", "SLC25A37", "PROK2", "BASP1", "CSF3R", # Neutrophils 10 markers (total until here 82)
                     "CORO1B", "CD27", "RTKN2", "ARID5B", "TRAC", "ITGB1", "TRBC2", "TRBC1", "IFI44L", "XAF1", "MX1", "IFIT3", "STAT1", "EIF2AK2", # CD4 T 14 markers (total until here 96)
                     "LST1", "FCGR3A", "AIF1", "MS4A7", "C1QA", "COTL1", "CFD", "PSAP", "FTL", # Non-classical Monocytes 9 markers (total until here 105)
                     "LTF", "LCN2", "CAMP", "RETN", "DEFA4", "CD24", "PGLYRP1", "MMP8", # Immature Neutrophils 8 markers (total until here 113)
                     "AC02O656.1", "TNFAIP2", "MPEG1", "CPVL", "FGL2", "CSTA", "CYBB", "MTRNR2L12", "CTSS", # Monocytes (ZFP36L2) 9 markers (total until here 122)
                     "STMN1", "MKI67", "TYMS", "PCNA", "TUBA1B", "TUBB", "DUT", "HMBG2", "HMGN2", "HIST1H4C", # Prol. T cells 10 markers (total until here 132)
                     "FCER1A", "HLA-DRB5", "HLA-BRB1", # mDCs 3 markers (total until here 135)
                     "JCHAIN", "IGHG1", "IGHG3", "MYB1", "IGHA1", "IGLC7", "IGLC2", "IGHM", "IGLC3", "IGKC", # Plasmablast 10 markers (total until here 145)
                     "PPBP", "PF4", "NRGN", "GNG11", "TUBB1", "CAVIN2", "HBG2", "MYL9", "GP9", "CLU", # Megakaryocytes 10 markers (total until here 155)
                     "HBA2", "HBA1", "HBD", "ALAS2", "HBM", "AHSP", "CA1", "SNCA", "HBB", # Erythroblasts 9 markers (total until here 164)
                     "SOX4", "PRSS57", "SPINK2", "AVP", "EGFL7", "CDK6", "ANKRD28", "SERPINB1", "AREG", # Undefined 9 markers (total until here 173)
                     "ITM2C", "PLD4", "SERPINF1", "TCF4", "LILRA4", "IRF8", "IRF7", "PTGDS", "PPP1R14B", # pDCs 9 markers (total until here 182)
                     "CMTM2", "SOD2", "PLEK" # CD8 3 markers (total until here 185)
)
marker.list <- list(
  c("CCR7", "LEHF1", "FHIT", "PIK3IP1", "LDHB", "MAL", "ADTRP", "EEF1B2", "RPS3A", "TCF7", "LTB", "IL7R", "AQP3", "JUNB", "TRADD", "RPSA", "RPLP0", "RPS6"),
  c("CORO1B", "CD27", "RTKN2", "ARID5B", "TRAC", "ITGB1", "TRBC2", "TRBC1", "IFI44L", "XAF1", "MX1", "IFIT3", "STAT1", "EIF2AK2"),
  c("VCAN", "LYZ", "S100A9", "S100A12", "LGALS1", "MAFB", "CD14", "PLBD1", "FCN1", "S100A8"),
  c("GZMH", "CCL5", "NKG7", "FGFBP2", "CD8A", "GZMA", "GNLY", "CST7", "CD8B", "CCL4", "GZMK", "DUSP2", "LYAR", "KLRG1", "IL32", "KLRB1", "CMTM2", "SOD2", "PLEK"),
  c("GZMB", "PRF1", "SPON2", "CLIC3", "CTSW", "HOPX"),
  c("IFII27", "IFITM3", "APOBEC3A", "TYMP", "ISG15", "CST3", "LGALS2", "IFI6"),
  c("CD79A", "MS4A1", "CD74", "CD79B", "HLA-DRA", "TCL1A", "HLA-DQB1", "HLA-DPA1", "HLA-DQA1", "HLA-DPB1"),
  c("S100B", "RPS5", "NELL2", "NOSIP"),
  c("FCGR3B", "NAMPT", "IFITM2", "CXCL8", "NEAT1", "G0S2", "SLC25A37", "PROK2", "BASP1", "CSF3R"),
  c("LST1", "FCGR3A", "AIF1", "MS4A7", "C1QA", "COTL1", "CFD", "PSAP", "FTL"),
  c("LTF", "LCN2", "CAMP", "RETN", "DEFA4", "CD24", "PGLYRP1", "MMP8"),
  c("AC02O656.1", "TNFAIP2", "MPEG1", "CPVL", "FGL2", "CSTA", "CYBB", "MTRNR2L12", "CTSS"),
  c("STMN1", "MKI67", "TYMS", "PCNA", "TUBA1B", "TUBB", "DUT", "HMBG2", "HMGN2", "HIST1H4C"),
  c("FCER1A", "HLA-DRB5", "HLA-BRB1"),
  c("JCHAIN", "IGHG1", "IGHG3", "MYB1", "IGHA1", "IGLC7", "IGLC2", "IGHM", "IGLC3", "IGKC"),
  c("PPBP", "PF4", "NRGN", "GNG11", "TUBB1", "CAVIN2", "HBG2", "MYL9", "GP9", "CLU"),
  c("HBA2", "HBA1", "HBD", "ALAS2", "HBM", "AHSP", "CA1", "SNCA", "HBB"),
  c("SOX4", "PRSS57", "SPINK2", "AVP", "EGFL7", "CDK6", "ANKRD28", "SERPINB1", "AREG"),
  c("ITM2C", "PLD4", "SERPINF1", "TCF4", "LILRA4", "IRF8", "IRF7", "PTGDS", "PPP1R14B"))

names(marker.list) <- c("CD4T.conventional1", "CD4T.conventional2", "Classical.Mono", "CD8T.EffectorMemory", "NK", "InflamMono", "B", "CD8T.N.NM", "Neutrophil",
                        "Non.Classical", "Immature.Neutro", "MonoZFP36L2", "Prol.T", "MDCs", "Plasmablasts", "Mega", "Erythroblast", "Undefined", "PDCs")

# Dotplots for the individual celltypes
for(i in 1:length(marker.list)){
  
  pdf(paste0(OUTPUT_DIR, "celltypes/dotplot_", names(marker.list)[i], ".pdf"), width = 10, height = 8)
  print(DotPlot(data,features=marker.list[[i]],cols="RdBu")+coord_flip() +
          theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
          ggtitle(label = paste0("PBMC Marker ", names(marker.list)[i])))
  dev.off()
}

# Dotplots of the marker genes more general
pdf(paste0(OUTPUT_DIR, "celltypes/dotplot_PBMC_markers.pdf"), width = 12, height = 9)
print(DotPlot(data,features=markers_pbmc,cols="RdBu")+coord_flip() +
        theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
        ggtitle(label = "PBMC Marker Genes"))
dev.off()

pdf(paste0(OUTPUT_DIR, "celltypes/dotplot_Neutrophil_markers.pdf"), width = 12, height = 9)
print(DotPlot(data,features=markers_neutrophil,cols="RdBu")+coord_flip() +
        theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
        ggtitle(label = "Neutrophil Marker Genes"))
dev.off()


# 2. Dotplots for top 20 markers per cellcluster

Idents(data) <- 'seurat_clusters'
markers <- FindAllMarkers(data, logfc_threshold = 0)

write.csv(markers, file = paste0(OUTPUT_DIR, 'celltypes/markers.csv'))

for (i in unique(data@meta.data$seurat_clusters)) {
  stopifnot(i %in% markers$cluster)
  DEG <- markers %>% filter(cluster == i) %>% arrange(p_val_adj) %>% head(20) %>% pull(gene)
  
  # Make dotplot
  pdf(paste0(OUTPUT_DIR, 'celltypes/dotplot_cluster_', i, '.pdf'), width = 10, height = 8)
  print(DotPlot(data, features = DEG, cols = "RdBu") +
          coord_flip() +
          theme(legend.position = "top", legend.title = element_blank()) +
          labs(x = "", y = "", title = 'markers')
  )
  dev.off()
}


# 3. Map cellcluster to celltypes based on plots from above
# Add the information in metadata and produce UMAP for the celltypes
data@meta.data$celltype <- recode(data@meta.data$seurat_clusters,
                                  `0` = 'Naive CD4+ T',
                                  `1` = 'NK',
                                  `2` = 'CD8+ T',
                                  `3` = 'B',
                                  `4` = 'Memory CD4+ T',
                                  `5` = 'Memory CD4+ T',
                                  `6` = 'CD14+ Monocyte',
                                  `7` = 'B',
                                  `8` = 'Unknown',
                                  `9` = 'Megakaryocyte',
                                  `10` = 'NK',
                                  `11` = 'CD16+ Monocytes',
                                  `12` = 'pDC',
                                  `13` = 'CD14+ Monocyte'
)

pdf(paste0(OUTPUT_DIR, "celltype_UMAP_largefont.pdf"), width = 18, height = 10)
print(DimPlot(data, group.by = "celltype", pt.size = 0.5, reduction = "umap", label = T, label.size = 6) +
        ggtitle(label = "UMAP of cell types") +
        NoLegend() +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

save.image(file = paste0(OUTPUT_DIR, 'data.RData'))

# 3. DE.R ---------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(openxlsx)

load('/vol/projects/mzoodsma/influenza/scRNAseq/analysis/output/data.RData')

# DEG between the conditions within each celltype
# Output is saved as .csv for each celltype and as .xlsx for all celtypes together

nr_DEGs <- matrix(0, nrow = length(unique(data@meta.data$celltype)), ncol = 2)
colnames(nr_DEGs) = c('Downregulated', 'Upregulated')
rownames(nr_DEGs) = unique(data@meta.data$celltype)
Idents(data) <- 'vaccination'

total <- list()

print(head(data@meta.data))
for (cell in unique(data@meta.data$celltype)) {
  
  sub <- subset(data, celltype == cell)
  
  # Use wilcoxon test for now. MAST when age and gender of the samples are known
  # Set ident.1 as after vaccination, and ident.2 as before vaccination. This way,
  # logFC values > 0 indicate upregulated genes AFTER vaccination.
  res <- FindMarkers(sub, ident.1 = 'after', ident.2 = 'before', logfc.threshold = 0)
  
  res$significance <- ifelse(res$p_val_adj < 0.05, TRUE, FALSE)
  res$direction <- ifelse(res$avg_log2FC > 0, 'Upregulated', 'Downregulated')
  res %>% arrange(p_val_adj)
  res$gene <- rownames(res)
  
  # Export
  write.csv(res, paste0('output/DEG/', cell, '_DEG.csv'))
  total[[cell]] <- res
  
  # Count number of sig genes and export this as well
  sig.res <- res %>% filter(p_val_adj < 0.05)
  nr_DEGs[cell, 'Downregulated'] <- sig.res %>% filter(avg_log2FC < 0) %>% nrow()
  nr_DEGs[cell, 'Upregulated'] <- sig.res %>% filter(avg_log2FC > 0) %>% nrow()
  
  # Volcano plot per celltype
  # Top X significant labeled
  pdf(paste0(OUTPUT_DIR, 'DEG/', cell, '_volcano.pdf'), width = 5, height = 5)
  print(ggplot() +
          geom_point(data = res %>% filter(significance == F) , aes(x = avg_log2FC, y = -log10(p_val_adj)), color = 'lightgray', show.legend = F) +
          geom_point(data = res %>% filter(significance == T) , aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
          labs(x = 'Log-fold change after / before vaccination', y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = paste0(cell)) +
          scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
          geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
          geom_vline(xintercept = 0.0, lty = 'dotted') +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggrepel::geom_label_repel(data = head(x = res, 15), aes(label = gene, x = avg_log2FC, y = -log10(p_val_adj)), box.padding = 0.5, show.legend = F),
        
  )
  dev.off()
  
}

write.csv(paste0(OUTPUT_DIR, 'DEG/nr_DEGs.csv'), x = nr_DEGs)
write.xlsx(file = paste0(OUTPUT_DIR, 'DEG/DEG_influenza.xlsx'), total)
# 4. cell_proportion.R ----------------------------------------
# Testing if different conditions have different cell proportions
# Author: Martin Grasshoff
# Date: 20.4.2020
# Result Object: None. All results are saved in csv files or pdf plots.
# Saved Objects:
# DirichletTestAdjusted.csv          Results of our regression
# DirichletTestAdjustedSig.csv       Results filtered down to only include statistically significant ones
# AllConditions_CellTypes.pdf        Boxplot of all conditions by cell type
# Condition_Critical_CellTypes.pdf   Boxplots for all critical samples by cell type. One image gets created per condition

rm(list = ls())
load('/vol/projects/mzoodsma/influenza/scRNAseq/analysis/output/data.RData')

TestingCellProportions <- function(data, identity.to.compare, condition.to.compare, output_path, type, celltypes, PlotSigIndicators = TRUE){
  
  # Libraries
  library(Seurat)
  library(ggplot2)
  library(ggsignif)
  library(reshape2)
  if(!require("DirichletReg")) install.packages("DirichletReg")
  
  # Reshaping the data for the Dirichlet-Test
  # Descriptive Statistics
  # Mean, Variance and Coefficient of Variance per Cell Type
  DescriptionCelltype <- function(SeuratObject, identity.to.compare, condition.to.compare, celltypes){
    
    # Result object
    result <- list()
    
    # Subsetting to the relevant subgroup
    data <- FetchData(SeuratObject, vars = identity.to.compare)
    data.subgroup <- SeuratObject[, which(x = data == condition.to.compare)]
    
    # Calculating the percentages per cell type per group
    result$CellsPerType <- table(as.character(data.subgroup$celltype))
    result$PercentagesPerType <- result$CellsPerType/dim(data.subgroup)[2]
    
    # Calculating the percentages per cell type per patient
    patients <- unique(data.subgroup$sample)
    patients <- patients[order(patients)]
    print(patients)
    CellsPerTypePerPatient <- c()
    PercentagesPerTypePerPatient <- c()
    values.for.MeanVariance <- c()
    for(person in patients){
      # Not all patients have all cell types at all time points.
      # With this dummy, we can use 0 for missing values
      result.dummy <- rep(0,length(celltypes))
      names(result.dummy) <- celltypes
      result.dummy <- result.dummy[sort(names(result.dummy))]
      
      # We subset to a sample ID (one person at one time point)
      data.person <- FetchData(data.subgroup, vars = "sample")
      data.person <- data.subgroup[, which(x = data.person == person)]
      
      # How many cells does a patient have of each type?
      result.preliminary <- table(as.character(data.person[[type]][,1]))
      result.preliminary <- result.preliminary[sort(names(result.preliminary))]
      
      result.dummy[names(result.dummy) %in% names(result.preliminary)] <- result.preliminary
      CellsPerTypePerPatient <- rbind(CellsPerTypePerPatient, result.dummy)
      rownames(CellsPerTypePerPatient)[dim(CellsPerTypePerPatient)[1]] <- paste0("Patient_", person)
      
      
      # Percentages
      result.dummy <- result.dummy/sum(result.dummy)
      PercentagesPerTypePerPatient <- rbind(PercentagesPerTypePerPatient, result.dummy)
      rownames(PercentagesPerTypePerPatient)[dim(PercentagesPerTypePerPatient)[1]] <- paste0("Patient_", person)
    }
    
    result$CellsPerPatient <- CellsPerTypePerPatient
    result$PercentagesPerPatient <- PercentagesPerTypePerPatient
    
    # Calculating the Mean, Variance and Coefficient of Variance per cell type
    result$MeanPerCelltype <- colMeans(PercentagesPerTypePerPatient)
    result$VarsPerCelltype <- apply(PercentagesPerTypePerPatient, 2, var)
    result$CoeffOfVariance <- result$VarsPerCelltype / result$MeanPerCelltype
    result$Condition <- condition.to.compare
    return(result)
  }
  
  CellTypePerCondition <- lapply(condition.to.compare, DescriptionCelltype,
                                 SeuratObject = data,
                                 identity.to.compare = identity.to.compare,
                                 celltypes = celltypes)
  
  
  # Boxplot to illustrate the variance
  data.boxplot.combined <- c()
  for(conditions in 1:length(condition.to.compare)){
    data.boxplot <- data.frame(patient = rownames(CellTypePerCondition[[conditions]]$PercentagesPerPatient), CellTypePerCondition[[conditions]]$PercentagesPerPatient, stringsAsFactors = FALSE)
    data.boxplot <- melt(data.boxplot, id.vars = "patient")
    
    pdf(paste0(output_path, "Condition_", condition.to.compare[conditions], "_CellTypes.pdf"))
    print(
      ggplot(data.boxplot, aes(x = variable, y = value)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    dev.off()
    
    # Combining the data to plot all conditions in one boxplot
    data.boxplot <- data.frame(data.boxplot, condition = condition.to.compare[conditions])
    data.boxplot.combined <- rbind(data.boxplot.combined, data.boxplot)
  }
  
  pdf(paste0(output_path, "AllConditions_CellTypes.pdf"), width = 10, height = 4)
  print(
    ggplot(data = data.boxplot.combined, aes(x = condition, y = value))+
      geom_boxplot(aes(fill = condition), outlier.shape = NA)+
      facet_grid(.~variable)+
      geom_jitter(shape=16, width = 0.2)+
      theme_classic()+
      xlab("cells")+ylab("cell proportion")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  
  # Performing the Dirichlet-Tests for each Cell Type
  # We loop over all conditions by prepending them with a LETTER
  
  results.all <- c() # We rbind our results to this variable
  for(i in 1:length(as.character(unique(data[[identity.to.compare]])[,1]))){
    # What are the LETTERS with which we prepend the conditions
    prepend <- LETTERS[unique(c(i:length(as.character(unique(data[[identity.to.compare]])[,1])), 1:i))]
    print(cat("What we prepend with: ", prepend))
    
    data.test <- data.frame(celltype  = as.character(data.boxplot.combined$variable),
                            sample    = as.character(data.boxplot.combined$patient),
                            value     = data.boxplot.combined$value,
                            other     = 1 - data.boxplot.combined$value,
                            condition = as.character(data.boxplot.combined$condition),
                            stringsAsFactors = FALSE)
    celltypes = unique(data.test$celltype)
    
    conditions <- unique(data.test$condition)
    conditions <- conditions[order(conditions)]
    
    for(prepending in 1:length(prepend)){
      data.test[data.test$condition == conditions[prepending], 5] <- paste0(prepend[prepending],conditions[prepending])
    }
    
    # Vector of all states in the identity to compare
    condition.vector <- unique(data.test$condition)
    condition.vector <- condition.vector[order(condition.vector)]
    comparisons <- paste0("condition", condition.vector)
    dat.test <- NULL
    for(type in celltypes){
      try({
        print(type)
        data.type <- data.test[data.test$celltype == type,]
        data.type <- rbind(subset(data.type, condition == condition.vector[1]),
                           subset(data.type, condition == condition.vector[2]),
                           subset(data.type, condition == condition.vector[3]),
                           subset(data.type, condition == condition.vector[4]))
        data.type$Smp <- DR_data(data.type[,3:4])
        res <- DirichReg(Smp~condition, data.type, model = "alternative", base = 1, verbosity = 0, control = list(iterlim = 5000))
        x <- summary(res)
        
        for(comparison in comparisons){
          if(comparison %in% rownames(x$coef.mat)){
            dat.test <- rbind(dat.test, c(type, x$coef.mat[comparison, 3:4], comparison))
            # Renaming the comparisons
            dat.test[dat.test[,4] == comparison,4] <- paste0(substr(condition.vector[1], start = 2, stop = nchar(condition.vector[1])), "Vs", substr(comparison, start = 11, stop = nchar(comparison)))
          }
        }
      })
    }
    
    dat.test <- as.data.frame(dat.test, stringsAsFactors = FALSE)
    colnames(dat.test) <- c("cell_type","Z","Pr","comparison")
    
    # Saving all results
    results.all <- rbind(results.all, dat.test)
  }
  # Renaming all the comparisons and removing duplicated ones
  old.comparisons <- results.all$comparison
  old.comparisons <- strsplit(old.comparisons, "Vs")
  old.comparisons <- lapply(old.comparisons, sort)
  new.comparisons <- lapply(old.comparisons, paste0, collapse = "Vs")
  new.comparisons <- unlist(new.comparisons)
  results.all$comparison <- paste0(results.all$cell_type, new.comparisons)
  results.all <- results.all[!duplicated(results.all$comparison),]
  
  # Bonferroni correcting the p-values
  results.all$sig <- "ns"
  results.all$PrAdjust <- p.adjust(results.all$Pr, method = "bonferroni")
  results.all$sig[results.all$PrAdjust < 0.05]="*"
  results.all$sig[results.all$PrAdjust < 0.01]="**"
  results.all$sig[results.all$PrAdjust < 0.001]="***"
  
  # Writing all results
  write.csv(results.all, paste0(output_path, "DirichletTestAdjusted.csv"), row.names = FALSE)
  # Saving only results that significant
  results.all.sig <- subset(results.all, sig != "ns")
  write.csv(results.all.sig, paste0(output_path, "DirichletTestAdjustedSig.csv"), row.names = FALSE)
}

# The data has to be read in before hand.

identity.to.compare <- 'vaccination'
condition.to.compare <- as.character(levels(as.factor(data@meta.data$vaccination)))
type <- 'celltype'
celltypes <- unique(as.character(data@meta.data$celltype))

output_path <- paste0(OUTPUT_DIR, 'cell_proportions/')
TestingCellProportions(data,
                       identity.to.compare,
                       condition.to.compare,
                       output_path,
                       type,
                       celltypes,
                       PlotSigIndicators = TRUE)
# 5. enrichment.R --------------------------------
rm(list = ls())

library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(AnnotationDbi)
library(ggpubr)
library(magrittr)
library(Seurat)

# celltypes <- c('Naive CD4+ T', 'Memory CD4+ T', 'CD8+ T', 'NK', 'CD14+ Monocyte', 'CD16+ Monocytes', 'B', 'pDC', 'Megakaryocyte', 'Unknown')
# Remove celltypes without DEG. These cause errors when trying to calculate enrichment
celltypes <- c('Naive CD4+ T', 'Memory CD4+ T', 'CD14+ Monocyte', 'NK', 'B')

#celltypes <- c('Memory CD4+ T')
for (celltype in celltypes){
  
  fname <- paste0('output/DEG/', celltype, '_DEG.csv')
  
  df <- read.csv(fname, header=T, row.names=1)
  up <- df %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0) %>% row.names()
  down <- df %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0) %>% row.names()
  
  print(paste0('Up: ', length(up)))
  print(paste0('Down: ', length(down)))
  
  up.GO = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)
  down.GO = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)
  
  enr.GO.up <- enrichGO(gene = up.GO,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  
  enr.GO.down <- enrichGO(gene = down.GO,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
  
  up.KEGG = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)
  down.KEGG = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)
  
  enr.KEGG.up <- enrichKEGG(gene = up.KEGG,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 1,
                            qvalueCutoff = 0.05)
  
  enr.KEGG.down <- enrichKEGG(gene = down.KEGG,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 1,
                              qvalueCutoff = 0.05)
  
  res = list('GO upreg' = enr.GO.up,
             'GO downreg' = enr.GO.down,
             'KEGG up' = enr.KEGG.up,
             'KEGG down' = enr.KEGG.down)
  
  openxlsx::write.xlsx(x = res, file = paste0('output/', celltype, '_enrichment_results.xlsx'))
  
  pdf(paste0('output/', celltype, 'enr.pdf'), width = 15, height = 15)
  print(merge_result(list('Upregulated GO' = enr.GO.up,
                          'Downregulated GO' = enr.GO.down,
                          'Upregulated KEGG' = enr.KEGG.up,
                          'Downregulated KEGG' = enr.KEGG.down)) %>%
          clusterProfiler::dotplot(., showCategory = 10) +
          theme(text = element_text(size =12),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
}

# 6. boxplot_specific_genes.R -------------------------------
rm(list = ls())

load('output/data.RData')

df <- FetchData(object = data, vars = c('FOS', 'JUN', 'NFKBIA', 'vaccination', 'celltype'), slot = 'data')

theme_set(theme_classic() + theme(text = element_text(size = 12)))

df$vaccination <- factor(df$vaccination, levels = c('before', 'after'))
# CD14+ Monocytes
sub <- df %>% filter(celltype == 'CD14+ Monocyte')
pdf('output/FOS_CD14_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = FOS)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'FOS (CD14+ Monocyte)', y = 'Normalised gene expression')
dev.off()

pdf('output/JUN_CD14_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = JUN)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'JUN (CD14+ Monocyte)', y = 'Normalised gene expression')
dev.off()

pdf('output/NFKBIA_CD14_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = NFKBIA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'NFKBIA (CD14+ Monocyte)', y = 'Normalised gene expression')
dev.off()


# CD16+ Monocytes
sub <- df %>% filter(celltype == 'CD16+ Monocytes')
pdf('output/FOS_CD16_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = FOS)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'FOS (CD16+ Monocyte)', y = 'Normalised gene expression')
dev.off()

pdf('output/JUN_CD16_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = JUN)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'JUN (CD16+ Monocyte)', y = 'Normalised gene expression')
dev.off()

pdf('output/NFKBIA_CD16_Mono.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = NFKBIA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'NFKBIA (CD16+ Monocyte)', y = 'Normalised gene expression')
dev.off()

# CD4+ Memory T cells
sub <- df %>% filter(celltype == 'Memory CD4+ T')
pdf('output/FOS_memory_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = FOS)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'FOS (Memory CD4+ T)', y = 'Normalised gene expression')
dev.off()

pdf('output/JUN_memory_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = JUN)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'JUN (Memory CD4+ T)', y = 'Normalised gene expression')
dev.off()

pdf('output/NFKBIA_memory_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = NFKBIA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'NFKBIA (Memory CD4+ T)', y = 'Normalised gene expression')
dev.off()

# Naive CD4+
sub <- df %>% filter(celltype == 'Naive CD4+ T')
pdf('output/FOS_naive_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = FOS)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'FOS (Naive CD4+ T)', y = 'Normalised gene expression')
dev.off()

pdf('output/JUN_naive_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = JUN)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'JUN (Naive CD4+ T)', y = 'Normalised gene expression')
dev.off()

pdf('output/NFKBIA_naive_T.pdf', width = 4, height = 4)
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = NFKBIA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'NFKBIA (Naive CD4+ T)', y = 'Normalised gene expression')
dev.off()
# 7. singeR.R --------------------------------------------

rm(list = ls())

library(Seurat)
library(SingleR)
library(celldex)
library(cowplot)
library(ggplot2)

ref.data1 <- HumanPrimaryCellAtlasData()
ref.data2 <- BlueprintEncodeData()
ref.data3 <- DatabaseImmuneCellExpressionData()
ref.data4 <- NovershternHematopoieticData()
ref.data5 <- MonacoImmuneData()

# Not included as ref.data:
#MouseRNAseqData()
#ImmGenData()

load('/vol/projects/mzoodsma/influenza/scRNAseq/analysis/output/data.RData')

annot1 <- SingleR(test = as.SingleCellExperiment(data), ref = ref.data1, labels = ref.data1$label.main, )
annot2 <- SingleR(test = as.SingleCellExperiment(data), ref = ref.data2, labels = ref.data2$label.main)
annot3 <- SingleR(test = as.SingleCellExperiment(data), ref = ref.data3, labels = ref.data3$label.main)
annot4 <- SingleR(test = as.SingleCellExperiment(data), ref = ref.data4, labels = ref.data4$label.main)
annot5 <- SingleR(test = as.SingleCellExperiment(data), ref = ref.data5, labels = ref.data5$label.main)

annot <- data.frame(HPCAD_labels = annot1$labels,
                    BED_labels   = annot2$labels,
                    DICED_labels = annot3$labels,
                    NHD_labels   = annot4$labels, 
                    MID_labels   = annot5$labels,
                    row.names = rownames(data@meta.data))

data <- AddMetaData(object = data, metadata = annot)

pdf(paste0(OUTPUT_DIR, 'singleR.pdf'), width = 30, height = 20)
plot_grid(
  DimPlot(object = data, group.by = 'HPCAD_labels', reduction = 'umap', pt.size = 0.5) + ggtitle('Human Primary Cell Atlas'),
  DimPlot(object = data, group.by = 'BED_labels', reduction = 'umap', pt.size = 0.5) + ggtitle('Blueprint Encode'),
  DimPlot(object = data, group.by = 'DICED_labels', reduction = 'umap', pt.size = 0.5) + ggtitle('Database of Immune cell expression'),
  DimPlot(object = data, group.by = 'NHD_labels', reduction = 'umap', pt.size = 0.5) + ggtitle('Hematopoietic database'),
  DimPlot(object = data, group.by = 'MID_labels', reduction = 'umap', pt.size = 0.5) + ggtitle('Immune cell populations')
)
dev.off()

# run in bash ----------------------------
#!/bin/bash

#$ -N run_R
#$ -l arch=linux-x64
#$ -b n
#$ -i /dev/null
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -q all.q


export RSCRIPT="/vol/biotools/F29/R-4.0.2/bin/Rscript"
$RSCRIPT ${1}

# run all in bash ---------------------------
#!/bin/bash

#$ -N run_R
#$ -l arch=linux-x64
#$ -b n
#$ -i /dev/null
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -q all.q

# Run all the R scripts to completely reproduce the results

export RSCRIPT="/vol/biotools/F29/R-4.0.2/bin/Rscript"
$RSCRIPT load_data.R
$RSCRIPT annotation.R
$RSCRIPT DE.R
$RSCRIPT cell_proportions.R
$RSCRIPT enrichment.R

