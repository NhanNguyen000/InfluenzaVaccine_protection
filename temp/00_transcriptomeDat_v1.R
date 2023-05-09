library(Seurat)
library(tidyverse)
library(patchwork)
# DATA from Martijn -------------------
# ATACseq data
/vol/projects/CIIM/Influenza/supplemental_data/Wimmers_et_al/tiv_scatac.h5ad.gz
#ATACseq from Wimmers et al. - CHECK THE CODE: /vol/projects/BIIM/MMR_vaccine/scATACseq/ArchR/Code,
# there are 2 main R package for ATACseq: Arch and SigNet (?), can ask Javi or Wenchao for ATACseq

# scRNAseq from Duesseldorf
/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData  
# marker.list is based on accumulated previous information (no specific papers)
# the marker table use differnt cutoff, - check the code at DE.R 

# iMED data (publish in Nature)- check the paper
/vol/projects/CIIM/Influenza/iMED/transcriptomic/
# for the iMED bulkRNA and scRNA (check the R file --> .rds files)

  
# Data from Samya --------------
Came across this influenza vaccination single cell atlas. Its only of 6 individuals, single cell, but they have quite a lot of cells i think. Some 112???662 cells. Maybe useful for looking for markers for minor cell populations???.
https://onlinelibrary.wiley.com/doi/10.1002/jmv.28174

https://onlinelibrary.wiley.com/doi/10.1002/jmv.28174

# # ATACseq data ------------
iMED scATACseq raw data: /vol/projects/CIIM/Influenza/iMED/single_cell/onlyATAC
iMED scATACseq demultiplexing results: /vol/projects/CIIM/Influenza/iMED/single_cell/onlyATAC/demultiplex/bwamem .
Pools that can be demultiplexed: FT4C (FT3C), FT4D (FT3D), FT6A, FT6C, FT6E

# INDIA vs dutch cohort. ----------------
/vol/projects/CIIM/INDIA/analysis/RNA/output/data.Rdata  
Different timepoints, stimulations for both populations. Metadata is all in there

# Dutch vs india ATACseq----------------
/vol/projects/CIIM/INDIA/analysis/ATAC/output/data.RData  Dutch vs india ATACseq
Baseline only, metadata in there also

# Associations of 300BCG and 500FG metabolites to cytokines and cellcounts in the subfolders here.
/vol/projects/CIIM/Influenza/iMED/supplemental_data  
# Associations of 300BCG and 500FG metabolites to cytokines and cellcounts in the subfolders here. Look in output/assoc_met_*.csv

dat <- read.csv("/vol/projects/CIIM/Influenza/iMED/supplemental_data/300BCG/output/assoc_met_cellcounts.csv")
dat <- read.csv("/vol/projects/CIIM/Influenza/iMED/supplemental_data/300BCG/output/assoc_met_cytokines.csv")
dat <- read.csv("/vol/projects/CIIM/Influenza/iMED/supplemental_data/300BCG/output/assoc_met_prot.csv")

C4H8O3_annot <- iMED$metabolite_annot %>% filter(Formula == "C4H8O3")

C4H8O3_assoc <- list()
assocTypes <- c("assoc_met_cellcounts.csv", "assoc_met_cytokines.csv", "assoc_met_prot.csv")
for (assocType in assocTypes) {
  dat <- read.csv(paste0("/vol/projects/CIIM/Influenza/iMED/supplemental_data/300BCG/output/", 
                         assocType))
  
  C4H8O3_assoc[["300BCG"]][[assocType]] <- dat %>% 
    filter(Top.annotation.ids %in% C4H8O3_annot$CompoundID) %>% 
    filter(p.value < 0.05)
}

assocTypes <- c("assoc_met_cellcount.csv", "assoc_met_cytokine.csv", "assoc_met_prot.csv")
for (assocType in assocTypes) {
  dat <- read.csv(paste0("/vol/projects/CIIM/Influenza/iMED/supplemental_data/500FG/output/", 
                         assocType)
  
  C4H8O3_assoc[["500FG"]][[assocType]] <- dat %>% 
    filter(Top.annotation.ids %in% C4H8O3_annot$CompoundID) %>% 
    filter(p.value < 0.05)
}

C4H8O3_annot <- iMED$metabolite_annot %>% filter(Formula == "C4H8O3")
b <- dat %>% filter(Top.annotation.ids %in% C4H8O3_annot$CompoundID) %>% 
  filter(p.value < 0.05)






#  scRNAseq from Duesseldorf ------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData")
a <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")
b<- read.csv("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/samplesheet_mkfastq.csv")

# plot of specific cell 

df <- FetchData(object = data, vars = c('FOS', 'JUN', 'NFKBIA', 'vaccination', 'celltype'), slot = 'data')
theme_set(theme_classic() + theme(text = element_text(size = 12)))
df$vaccination <- factor(df$vaccination, levels = c('before', 'after'))
sub <- df %>% filter(celltype == 'CD14+ Monocyte') # 'CD16+ Monocytes', 'Memory CD4+ T'
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = FOS)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'FOS (CD14+ Monocyte)', y = 'Normalised gene expression')

# test code ----
df <- FetchData(object = data, vars = c('CD83', 'vaccination', 'celltype'), slot = 'data')
df$vaccination <- factor(df$vaccination, levels = c('before', 'after'))
unique(df$celltype)
sub <- df %>% filter(celltype == 'pDC')
ggplot(sub) +
  geom_boxplot(aes(x = vaccination, y = CD83)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'CD83 (pDC)', y = 'Normalised gene expression')


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
  # write.csv(res, paste0('output/DEG/', cell, '_DEG.csv'))
  total[[cell]] <- res
  
  # Count number of sig genes and export this as well
  sig.res <- res %>% filter(p_val_adj < 0.05)
  nr_DEGs[cell, 'Downregulated'] <- sig.res %>% filter(avg_log2FC < 0) %>% nrow()
  nr_DEGs[cell, 'Upregulated'] <- sig.res %>% filter(avg_log2FC > 0) %>% nrow()
  
  # Volcano plot per celltype
  # Top X significant labeled
  # pdf(paste0(OUTPUT_DIR, 'DEG/', cell, '_volcano.pdf'), width = 5, height = 5)
  # print(ggplot() +
  #         geom_point(data = res %>% filter(significance == F) , aes(x = avg_log2FC, y = -log10(p_val_adj)), color = 'lightgray', show.legend = F) +
  #         geom_point(data = res %>% filter(significance == T) , aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
  #         labs(x = 'Log-fold change after / before vaccination', y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = paste0(cell)) +
  #         scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
  #         geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
  #         geom_vline(xintercept = 0.0, lty = 'dotted') +
  #         theme_classic() +
  #         theme(plot.title = element_text(hjust = 0.5)) +
  #         ggrepel::geom_label_repel(data = head(x = res, 15), aes(label = gene, x = avg_log2FC, y = -log10(p_val_adj)), box.padding = 0.5, show.legend = F),
  #       
  # )
  # dev.off()
  
}

# write.csv(paste0(OUTPUT_DIR, 'DEG/nr_DEGs.csv'), x = nr_DEGs)
# write.xlsx(file = paste0(OUTPUT_DIR, 'DEG/DEG_influenza.xlsx'), total)

# Testing if different conditions have different cell proportions
# Author: Martin Grasshoff; Date: 20.4.2020
# Result Object: None. All results are saved in csv files or pdf plots.
# Saved Objects:
# DirichletTestAdjusted.csv          Results of our regression
# DirichletTestAdjustedSig.csv       Results filtered down to only include statistically significant ones
# AllConditions_CellTypes.pdf        Boxplot of all conditions by cell type
# Condition_Critical_CellTypes.pdf   Boxplots for all critical samples by cell type. One image gets created per condition

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

## Nhan code -----------------------
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
DimPlot(data, reduction = "pca")
DimPlot(data, reduction = "umap")
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

a <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
identical(a, markers) # different?
Ident(data)
names(marker.list) # no dendric cell?

data2 <- RenameIdents(data2, marker.list)

# Conventional dendritic cell markers, based on https://hbctraining.github.io/scRNA-seq/lessons/08_SC_clustering_quality_control.html
FeaturePlot(data, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Plasmacytoid dendritic cell markers
FeaturePlot(data, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#  scRNAseq from iMED data ------------
list.files("/vol/projects/CIIM/Influenza/iMED/transcriptomic/")
data <- readRDS("/vol/projects/CIIM/Influenza/iMED/transcriptomic/iMED_published_scRNAseq.rds")
dat <- a

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(Seurat)
library(rstatix)
library(ggpubr)

data <- readRDS('/vol/projects/CIIM/Influenza/iMED/transcriptomic/scRNAseq/influenzaData_reanalysed.rds')
DimPlot(data, group.by = 'celltype') + scale_color_npg() + theme(aspect.ratio = 1)

# Metabolome genesets enriched at T3 vs T1 in TR
load('/vol/projects/CIIM/Influenza/iMED/transcriptomic/TR_T3vsT1ANDT4vsT1_enrichment.RData')

doAddModule <- function(seur, enr, padj.thr = 0.1) {
  # Add a module scores for all leading edge genes in the enr object to seur
  
  enr %<>% filter(padj <= padj.thr)
  l <- list()
  
  for (i in 1:nrow(enr)) {
    row <- enr[i, ]
    l[[row$pathway]] <- unlist(row$leadingEdge)
    
  }
  
  return(AddModuleScore(object = seur, ctrl = 5, features = l, name = paste0('enrscore_', enr$pathway)))
}

data <- doAddModule(data, mylist$TR_T3vsT1_metabolismEnrichment)

# Statistics between the celltypes for each of the pathways sep

metadata <- data@meta.data
pathway <- 'enrscore_REACTOME_METABOLISM_OF_RNA3'

doMake <- function(metadata, pathway, filename) {
  
  df <- metadata %>% 
    select(samptype, days, celltype, pathway) %>% 
    melt()
  
  # png(filename, width = 10, height = 5, res=  600, units = 'in')
  print(ggplot(df) + 
          geom_violin(aes(x = celltype, fill = celltype, y = value)) + 
          geom_boxplot(aes(x = celltype, fill = celltype, y = value), position = position_dodge(width = .9), width = .2, show.legend = F) +
          theme_classic() +
          theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
                plot.title = element_text(hjust = .5),
                legend.position = 'none') +
          labs(x = 'Celltype', y = 'Module score', title = pathway) +
          facet_grid('days') +
          scale_fill_npg() )
  
  
  # FeaturePlot(data, 
  #             features = pathway) + theme(aspect.ratio = 1) )
  dev.off()
  
}
mylist$TR_T3vsT1_metabolismEnrichment %>% filter(padj <0.1) %>% pull(pathway)


doMake(data@meta.data, 'enrscore_REACTOME_INTEGRATION_OF_ENERGY_METABOLISM1', 'scRNAseq/output/integration_energy_met.png')
doMake(data@meta.data, 'enrscore_REACTOME_METABOLISM_OF_CARBOHYDRATES2', 'scRNAseq/output/carbohydrate_metabolism.png')
doMake(data@meta.data, 'enrscore_REACTOME_METABOLISM_OF_RNA3', 'scRNAseq/output/rna_met.png')
doMake(data@meta.data, 'enrscore_REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE4', 'scRNAseq/output/tca_cycle.png')

# Nhan code --------------
data <- readRDS('/vol/projects/CIIM/Influenza/iMED/transcriptomic/scRNAseq/influenzaData_reanalysed.rds')
DimPlot(data, group.by = 'celltype') + scale_color_npg() + theme(aspect.ratio = 1)

unique(data$celltype)
df <- FetchData(object = data, vars = c('CD83', 'celltype'), slot = 'data')
ggplot(df) +
  geom_boxplot(aes(x = celltype, y = CD83)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'CD83 (across cell types)', y = 'Normalised gene expression') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

names(data@meta.data)
df <- FetchData(object = data, 
                vars = c('CD83', 'celltype',
                         "sampid", "samptype", "days", "sample", "condition", "cond"), 
                slot = 'data')

ggplot(df) +
  geom_boxplot(aes(x = celltype, y = CD83, fill = samptype)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'CD83 (across cell types)', y = 'Normalised gene expression') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# metadata from reclassification --------------------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

metadat_v2 <- metadat %>% filter(cohort == "iMED") %>%
  mutate(sampid = gsub("-", "", patientID)) %>% filter(sampid %in% df$sampid)

df2 <- df %>% full_join(metadat_v2 %>% select(sampid, matches("reclassify")))

ggplot(df2) +
  geom_boxplot(aes(x = celltype, y = CD83, fill = H1N1_reclassify)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'CD83 (across cell types)', y = 'Normalised gene expression') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
