rm(list = ls())
library(tidyverse)
library(caret)

get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[, groupType])) +
    geom_point(alpha=.5) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
    ggsci::scale_color_nejm() +
    labs(x = paste0("PC1 [",
                    round(100*summary(pca)$importance["Proportion of Variance", "PC1"],
                          digits = 2),
                    "%]"),
         y = paste0("PC2 [",
                    round(100*summary(pca)$importance["Proportion of Variance", "PC2"],
                          digits = 2),
                    "%]"),
         color = groupType) #+ stat_ellipse()
  
}

# load data --------------------------
load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# protein data -------------------------
proteinDat_metadata <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% select(group, name) %>% 
  inner_join(metadata_healthy)

proteinDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% filter(name %in% proteinDat_metadata$name) %>%
  select(-group) %>% column_to_rownames("name")

protein_pca <- predict(preProcess(x=proteinDat, method = "knnImpute"), # imput the NA value
                       proteinDat) %>% prcomp()
summary(protein_pca)$importance[1:3, 1:5]

# PCA plot with function
get.pca_plot(pca = protein_pca, 
             metadat = proteinDat_metadata, 
             groupType = "season")

get.pca_plot(pca = protein_pca, 
             metadat = proteinDat_metadata, 
             groupType = "responder")
# manual PCA plot
ggplot(data.frame(protein_pca$x), aes(PC1, PC2, color = proteinDat_metadata$season)) +
  geom_point(alpha=.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  ggsci::scale_color_nejm()

# metabolite data -------------------------
metaboliteDat_metadata <- mebo_Dat %>%
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% select(group, name) %>%
  inner_join(metadata_healthy)

metaboliteDat <- mebo_Dat %>%
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% filter(name %in% metaboliteDat_metadata$name) %>%
  select(-group) %>% column_to_rownames("name")

metabolite_pca <- metaboliteDat %>% select(-matches("_2")) %>% # with metabolites has double ionIdxs in iMED, remove the measures of the 2nd ionIdx
  prcomp()
summary(metabolite_pca)$importance[1:3, 1:5]

# PCA plot with function
get.pca_plot(pca = metabolite_pca, 
             metadat = metaboliteDat_metadata, 
             groupType = "season")
