library(Seurat)
library(tidyverse)
library(patchwork)
#  scRNAseq from Duesseldorf ------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData")
# Note: marker.list is based on accumulated previous information (no specific papers)
metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

#VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# subset only scRNAseq before vaccination
unique(data$vaccination)
dat2 <- subset(x = data, subset = vaccination == "before")

# DimPlot(dat2, reduction = "pca")
# DimPlot(dat2, reduction = "umap")
# DimPlot(dat2, group.by = 'celltype')
DimPlot(dat2, group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()

# plot of specific cell 
# FeaturePlot(dat2, reduction = "umap", 
#             features = c("CD83"), sort.cell = TRUE,
#             min.cutoff = 'q10', label = TRUE)
FeaturePlot(dat2, reduction = "umap", 
            features = c("CD83"), sort.cell = TRUE,
            min.cutoff = '1', label = TRUE)

# # boxplot
# df <- FetchData(object = data, vars = c('CD83', "sample", "vaccination", 'celltype'), slot = 'data') %>%
#   filter(vaccination ==  "before")
# 
# ggplot(df) +
#   geom_boxplot(aes(x = celltype, y = CD83)) + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 30, hjust = 1), 
#         plot.title = element_text(hjust = 0.5)) + 
#   labs(title = 'CD83', y = 'Normalised gene expression')
protein <- "CD83"
FeaturePlot(dat2, reduction = "umap", 
            features = protein, sort.cell = TRUE,
            min.cutoff = '1', label = TRUE)

RidgePlot(dat2, features = protein, group.by = "celltype") + NoLegend()
DotPlot(dat2, features = protein, group.by = "celltype") + RotatedAxis()
VlnPlot(dat2, features =protein, group.by = "celltype") + NoLegend()

protein <- "CD83"
FeaturePlot(dat2, reduction = "umap", 
            features = protein, sort.cell = TRUE,
            min.cutoff = '1', label = TRUE)
protein <- "CD83"
protein <- "CD79B"
protein <- "CD22"
protein <- "FCRL2"
protein <- "SIGLEC10"
protein <- "PARP1"
cowplot::plot_grid(
  VlnPlot(dat2, features =protein, group.by = "celltype") + NoLegend(),
  DotPlot(dat2, features = protein, group.by = "celltype") + RotatedAxis()
)



# iMED data (publish in Nature)- ----------------
data <- readRDS("/vol/projects/CIIM/Influenza/iMED/transcriptomic/iMED_published_scRNAseq.rds")

# subset only scRNAseq before vaccination
unique(data$days)
dat2 <- subset(x = data, subset = (days == "D0"))

length(unique(dat2$sample))
unique(dat2$sample)
unique(dat2$samptype)

# metadata for reclassification
# metadata from reclassification
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

VlnPlot(dat2, features = "CD83", group.by = "celltype") + NoLegend()
# get specific cell types
dat3 <- subset(x = dat2, subset =(celltype == "B" | celltype == "Plasmablast" | celltype == "mix"))

View(dat3@meta.data)
table(Idents(dat3))

a <- dat3@meta.data %>% mutate(patientID = gsub("I", "I-", sampid)) %>% 
  left_join(metadat %>% select(patientID, matches("reclassify")))

a %>% count(H1N1_reclassify, sample)
a %>% count(H3N2_reclassify, sample)
a %>% count(H1N1_reclassify, sample)

dat3@meta.data <- dat3@meta.data %>% mutate(patientID = gsub("I", "I-", sampid)) %>% 
  left_join(metadat %>% select(patientID, matches("reclassify")))

VlnPlot(dat3, features = "CD83", group.by = "celltype") + NoLegend()
VlnPlot(dat3, features = "CD83", group.by = "celltype", split.by = "H1N1_reclassify")
VlnPlot(dat3, features = "CD83", group.by = "celltype", split.by = "H3N2_reclassify")
VlnPlot(dat3, features = "CD83", group.by = "celltype", split.by = "B_reclassify")

# boxplot
df <- FetchData(object = dat2, 
                vars = c('CD83', "sample", "sampid", "samptype", "condition", 'celltype'), slot = 'data') 

plotDat <- df %>% mutate(patientID = gsub("I", "I-", sampid)) %>% 
  left_join(metadat %>% select(patientID, matches("reclassify"))) %>%
  filter(celltype %in% c("B", "Plasmablast", "mix"))

# compare re-classfiy groups
my_comparisons_v2.1 <- list( c("LL", "LH"), c("LL", "HL"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))

strain <- "H1N1_reclassify"
strain <- "H3N2_reclassify"
strain <- "B_reclassify"
cellType <- "B"
protein <- "CD83"


plotDat %>% 
  ggboxplot(x = strain, y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = my_comparisons_v2.1, method = "t.test")


plotDat %>% 
  ggboxplot(x = strain, y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

plotDat %>% 
  mutate(samptype = factor(samptype, levels = c("NR", "R"))) %>%
  ggboxplot(x = "samptype", y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = list( c("NR", "R")), method = "t.test")

# violin  plot
plotDat %>% 
  ggboxplot(x = strain, y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = my_comparisons_v2.1, method = "t.test")


plotDat %>% 
  ggviolin(x = strain, y = protein,
           paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

plotDat %>% 
  mutate(samptype = factor(samptype, levels = c("NR", "R"))) %>%
  ggviolin(x = "samptype", y = protein,
           paletter = "jco", add = "jitter") + facet_wrap(~celltype)+
  stat_compare_means(comparisons = list( c("NR", "R")), method = "t.test")

