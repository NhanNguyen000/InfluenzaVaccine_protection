library(readxl)
library(tidyverse)
library(openxlsx)

# functions ------------------------------
get.celltype_metadata <- function(rawDat) {
  metadata <- rawDat[1, -c(1:2)] %>% 
    t() %>% as.data.frame() %>%
    rownames_to_column(var="cell_type") %>%
    rename("subcell_type" = "V1")
  
  metadata$cell_type[grep("[...]", metadata$cell_type)] <- NA 
  metadata <- metadata %>% fill(cell_type)
  
  return(metadata)
}

get.celltype_data <- function(rawDat) {
  celltype_data <- rawDat
  colnames(celltype_data) <- c( "condition", celltype_data[1, -1])
  
  celltype_data <- celltype_data %>% 
    slice(-1) %>% drop_na("Donor ID") %>%
    fill(condition) %>% as.data.frame()
  
  celltype_data[, 3:53] <- lapply(celltype_data[, 3:53], as.numeric)
  
  return(celltype_data)
}
#  available samples ------------------------------
sheet_names <- c("Baseline ex vivo", "Baseline re-stimulated", 
                 "Visit 2 ex vivo", "Visit 2 re-stimulated")

flowCyto_metadata <- list()
flowCyto_data <- list()

for (sheet_name in sheet_names) {
  dat_temp <- read_xlsx("../flowCytometry/20230105_ZirFlu2019-2020_FlowCytometry_Valerie&Nhan.xlsx",
                        sheet = sheet_name)
  
  flowCyto_metadata[["2019"]][[sheet_name]] <- get.celltype_metadata(dat_temp)
  flowCyto_data[["2019"]][[sheet_name]] <- get.celltype_data(dat_temp) %>% 
    mutate(condition2 = factor(ifelse(condition == "Healthy", "healthy", "cirrhosis")))
  
  dat_temp <- read_xlsx("../flowCytometry/20230105_ZirFlu2020-2021_FlowCytometry_Valerie&Nhan.xlsx",
                        sheet = sheet_name)
  
  flowCyto_metadata[["2020"]][[sheet_name]] <- get.celltype_metadata(dat_temp)
  flowCyto_data[["2020"]][[sheet_name]] <- get.celltype_data(dat_temp) %>% 
    mutate(condition2 = factor(ifelse(condition == "Healthy", "healthy", "cirrhosis")))
}

library(ggplot2)
library(ggpubr)

get.ggboxplot <- function(dat, x_col, y_col) {
  ggboxplot(dat, x = x_col, y = y_col,
            color = "condition2", 
            paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test")
}

flowCyto_metadata <- flowCyto_metadata$`2019`
flowCyto_data <- flowCyto_data$`2019`
plot_list <- list()
for(celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  plot_list[[celltype]] <- cowplot::plot_grid(
    get.ggboxplot(flowCyto_data$`Baseline ex vivo`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Baseline re-stimulated`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Visit 2 ex vivo`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Visit 2 re-stimulated`,
                  x_col = "condition2", y_col = celltype),
    nrow = 1
  )
}

# pdf(file = "temp/20230206_cytoflowmetry_boxplot_2019.pdf", 
#     width = 12, height = 4 , onefile = TRUE)
# for (celltype in names(plot_list)) {
#   plot(plot_list[[celltype]])
# }
# dev.off()

plot_list$Th17
plot_list$Th1
plot_list$Th2
plot_list$`CD4+ Central Memory`
plot_list$`CD4+ Effector`
plot_list$`CD4+ IL-10 Tregs`
plot_list$`CD8+ Effector Memory`
plot_list$`CD8+ Effector`
plot_list$`CD8+ CD28+`
plot_list$`CD8+ CD38+`
plot_list$`CD4+ IL-10+`
plot_list$`CD8+ IFNy+`
plot_list$`CD8+ IL-10+`
plot_list$`B cells`
plot_list$`CD80 activated`
plot_list$`CD86 activated`
plot_list$Naive
plot_list$Memory
plot_list$`IgA switched`
plot_list$`IgG switched`
plot_list$`IgM non-switched`
plot_list$Bregs
plot_list$`IL35 Bregs`
plot_list$Transitional
plot_list$`Plasma cells`
plot_list$Plasmablasts

# check the expression based on re-classification -----------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis")))

flowCyto_data_total <- flowCyto_data %>% 
  lapply(function(x) x %>% 
           imap(~.x %>% mutate(measure = paste(.y))) %>%
           purrr::reduce(full_join)) %>%
  imap(~.x %>% mutate(season = paste(.y))) %>% 
  purrr::reduce(full_join)

plotDat <- flowCyto_data_total %>% 
  left_join(metadat, by = c("season", c("Donor ID" = "patientID")))

library(ggpubr)
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("HL", "HH"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"), c("LL", "LH", "HL", "HH"))


cellType <- "Th17"
# per measure
plotDat_temp <- plotDat %>% filter(measure = "Baseline ex vivo")
ggboxplot(plotDat_temp, x = "category", y = cellType,
          paletter = "jco", add = "jitter") + facet_wrap(~group) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# all measure
ggboxplot(plotDat, x = "category", y = cellType,
          paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(plotDat, x = "H1N1_reclassify", y = cellType,
          paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(plotDat, x = "H3N2_reclassify", y = cellType,
          paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(plotDat, x = "Bvictoria_reclassify", y = cellType,
          paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(plotDat, x = "Byamagata_reclassify", y = cellType,
          paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

# because ZirFlu have enough LL group in H1N1 
cellType <- "Th17"
plot_list <- list()
for (cellType in flowCyto_metadata$`2019`$`Baseline ex vivo`$subcell_type) {
  plot_list[[cellType]] <- cowplot::plot_grid(
    ggboxplot(plotDat, x = "category", y = cellType,
              paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
      stat_compare_means(comparisons = my_comparisons, method = "t.test"),
    
    ggboxplot(plotDat, x = "H1N1_reclassify", y = cellType,
              paletter = "jco", add = "jitter") + facet_grid(group~measure) + 
      stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
    
    ncol = 2
  )
}

pdf(file = "20230418_cellpropotion_ZirFlu_2seasons.pdf", width = 15, onefile = TRUE)
for (i in names(plot_list)) {
  plot(plot_list[[i]])
}
dev.off()

plot_list$Th17
plot_list$`B cells`
### check the association of CD83 to cell proportion ------------
# CD*3 - OID20565
plotDat_protein <- protein_normDat$ZirFlu %>% 
  select(OID20565) %>% rename("CD83" = "OID20565") %>%
  rownames_to_column("name") %>% 
  right_join(plotDat)

plotDat_protein %>%
  ggplot(aes(x = CD83, y = 'B cell')) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 7) +
  theme_classic() + ggtitle("HAI titer of the H1N1 strain")

plotDat_protein %>% filter(condition2 == "healthy") %>% 
  ggscatter(x = "CD83", y = "B cells", add = "reg.line", cor.coef = TRUE)
plotDat_protein %>% filter(condition2 == "healthy") %>% 
  ggscatter(x = "CD83", y = "Memory", add = "reg.line", cor.coef = TRUE)
plotDat_protein %>% filter(condition2 == "healthy") %>% 
  ggscatter(x = "CD83", y = "Naive", add = "reg.line", cor.coef = TRUE)
plotDat_protein %>% filter(condition2 == "healthy") %>% 
  ggscatter(x = "CD83", y = "IgA switched", add = "reg.line", cor.coef = TRUE)

# check correlation 
plotDat_CD83_healthy <- plotDat_protein %>% filter(condition2 == "healthy")

outcome_temp <- c()
for (cell in flowCyto_metadata$`2019`$`Baseline ex vivo`$subcell_type) {
  dat_temp <- plotDat_CD83_healthy %>% select(CD83, any_of(cell))
  cor_res <- cor_test(dat_temp)
  outcome_temp <- rbind(outcome_temp, cor_res)
}

outcome_temp %>% filter(p < 0.05)
cowplot::plot_grid(
  plotDat_protein %>% filter(condition2 == "healthy") %>% 
    ggscatter(x = "CD83", y = "B cells", add = "reg.line", cor.coef = TRUE),
  plotDat_protein %>% filter(condition2 == "healthy") %>% 
    ggscatter(x = "CD83", y = "CD21low", add = "reg.line", cor.coef = TRUE),
  plotDat_protein %>% filter(condition2 == "healthy") %>% 
    ggscatter(x = "CD83", y = "CD80 activated", add = "reg.line", cor.coef = TRUE),
  plotDat_protein %>% filter(condition2 == "healthy") %>% 
    ggscatter(x = "CD83", y = "Naive", add = "reg.line", cor.coef = TRUE)
)
plotDat_protein %>% filter(condition2 == "healthy") %>% 
  ggscatter(x = "CD83", y = "Naive", add = "reg.line", cor.coef = TRUE)

# check if protein & cell type associated -----------------
# get sample at baseline
sample_T1 <- ZirFlu$donorSamples %>% filter(time == "T1" & season == "2019")

# proteins
selected_proteins <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("TNFSF10", "CXCL8", "IL6", "IL17C", "IL17D", "TNFSF12"))
selected_proteins <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("CD83"))

protein_dat <- ZirFlu$protein_dat %>% select(selected_proteins$OlinkID) %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID) %>% 
  left_join(ZirFlu$donorSamples) %>% select(-c(probenID, season, time)) %>%
  relocate(patientID)

View(flowCyto_data$`Baseline ex vivo`)

# check the association - linear (check the code later)--------------------
library(rstatix)

# proteins
protein_dat <- ZirFlu$protein_dat %>% select(selected_proteins$OlinkID) %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID) %>%
  left_join(ZirFlu$donorSamples %>% 
              filter(season == "2019") %>% select(probenID, patientID))

outcome <- list()
for (protein in selected_proteins$OlinkID) {
  outcome_temp <- c()
  
  for (cell in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
    dat_temp <- protein_dat %>% select(patientID, any_of(protein)) %>%
      full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(cell))),
                by = c("patientID" = "Donor ID")) %>%
      column_to_rownames("patientID")
    
    cor_res <- cor_test(dat_temp)
    outcome_temp <- rbind(outcome_temp, cor_res)
  }
  outcome[[protein]] <- outcome_temp
}

# need to check all cells - no Th1, Th17, B cell show up
outcome_pval <- outcome %>% lapply(function(x) x %>% filter(p < 0.05))

plotList_1 <- list()
plotList_2 <- list()
for (celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  for (protein in selected_proteins$OlinkID) {
    plot_dat <- protein_dat %>% select(patientID, any_of(protein)) %>%
      full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(celltype))),
                by = c("patientID" = "Donor ID")) %>%
      left_join(ZirFlu$donorInfo %>% filter(season == "2019"))
    
    plotList_1[[celltype]][[protein]] <- plot_dat %>% 
      ggscatter(x = celltype, y = protein, 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson") +
      xlab(paste0(celltype, " percentage")) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
    
    plotList_2[[celltype]][[protein]] <- plot_dat %>% 
      ggscatter(x = celltype, y = protein, add = "reg.line") +
      stat_cor(method = "pearson") +
      xlab(paste0(celltype, " percentage")) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
    
  }
}

pdf(file = "temp/20230206_proteinT1_associCellType_T1.pdf", width = 15, onefile = TRUE)
for (celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_1[[celltype]]$OID20477, 
                       plotList_1[[celltype]]$OID20481,
                       plotList_1[[celltype]]$OID20563,
                       plotList_1[[celltype]]$OID20611,
                       plotList_1[[celltype]]$OID20631,
                       nrow = 1),
    cowplot::plot_grid(plotList_2[[celltype]]$OID20477, 
                       plotList_2[[celltype]]$OID20481,
                       plotList_2[[celltype]]$OID20563,
                       plotList_2[[celltype]]$OID20611,
                       plotList_2[[celltype]]$OID20631,
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()


plot_dat <- protein_dat %>% select(patientID, any_of(protein)) %>%
  full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(celltype))),
            by = c("patientID" = "Donor ID")) %>%
  left_join(ZirFlu$donorInfo %>% filter(season == "2019"))

plot_dat %>% 
  ggscatter(x = celltype, y = protein, 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") +
  xlab(paste0(celltype, " percentage")) + 
  ylab(paste0("log2_expression_", 
              ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))

protein <- "CD83"
celltype <- "B cells"
plot_dat %>% 
  ggscatter(x = celltype, y = protein, add = "reg.line") +
  stat_cor(method = "pearson") +
  xlab(paste0(celltype, " percentage")) + 
  ylab(paste0("log2_expression_", 
              ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
