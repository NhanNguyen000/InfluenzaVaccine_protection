rm(list = ls())
library(tidyverse)
library(caret)

# load data --------------------------
load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

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

# manual PCA plot-------------------------------------
pca_plotDat <- metabolite_pca$x %>% as.data.frame %>% rownames_to_column("name") %>% full_join(metaboliteDat_metadata)

## season ----------------------------------------
pcaSeason <- pca_plotDat %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = season), size = 3, alpha=.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  ggsci::scale_color_nejm() +
  labs(
    x = paste0("PC1 [", round(100*summary(metabolite_pca)$importance["Proportion of Variance", "PC1"], digits = 2), "%]"),
    y = paste0("PC2 [", round(100*summary(metabolite_pca)$importance["Proportion of Variance", "PC2"], digits = 2), "%]")) + 
  theme(legend.position = "top", text = element_text(size = 16))

# save the plot 
png("output/pcaMetabolite_perSeason.png")
pcaSeason
dev.off()


## season and sex --------------------
pca_plotDat %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = season, shape = sex), size = 2, alpha=.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  ggsci::scale_color_nejm() +
  labs(
    x = paste0("PC1 [", round(100*summary(metabolite_pca)$importance["Proportion of Variance", "PC1"], digits = 2), "%]"),
    y = paste0("PC2 [", round(100*summary(metabolite_pca)$importance["Proportion of Variance", "PC2"], digits = 2), "%]")) + 
  theme(legend.position = "top")

