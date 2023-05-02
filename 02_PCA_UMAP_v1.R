# check UMAP (not work) & PCA plot (could not separate categories) ------------------------------------
library(umap)
imed.umap <- umap(proteinDat_impute$iMED)
plot.iris(imed.umap)

metadata <- ZirFlu$donorSamples %>% full_join(ZirFlu$donorInfo) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  mutate(age_class = ifelse(age >60, ">60", "young"))

protein_pca <- predict(preProcess(x=ZirFlu$protein_dat, method = "knnImpute"),
                       ZirFlu$protein_dat) %>% prcomp()
summary(protein_pca)$importance[1:3, 1:5]
get.pca_plot(protein_pca, metadata, "condition")

metabolite_pca <- ZirFlu$metabolite_dat %>% prcomp()
summary(metabolite_pca)$importance[1:3, 1:5]
get.cowplot_pca(pca_dat = metabolite_pca, metadata)

get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[, groupType])) +
    geom_point(alpha=.5) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
    ggsci::scale_color_nejm()# +
  # labs(x = paste0("PC1 [", 
  #                 round(100*summary(pca)$importance["Proportion of Variance", "PC1"], 
  #                       digits = 2), 
  #                 "%]"),
  #      y = paste0("PC2 [", 
  #                 round(100*summary(pca)$importance["Proportion of Variance", "PC2"], 
  #                       digits = 2), 
  #                 "%]"), 
  #      color = groupType) +
  # stat_ellipse()
  # 
}

metadat_iMED_temp <- cohorts_dat$donorSample_all %>% 
  filter(name %in% rownames(proteinDat_impute$iMED)) %>% 
  left_join(cohorts_dat$HAItiter_all)
protein_pca <- proteinDat_impute$iMED %>% prcomp()
groupType = "category"

get.pca_plot(pca = protein_pca,
             metadat = metadat_iMED_temp, 
             groupType = "category")

ggplot(data.frame(protein_pca$x), aes(PC1, PC2, color = metadat_iMED_temp$category)) +
  geom_point(alpha=.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  ggsci::scale_color_nejm()


protein_pca_inflam <- proteinDat_impute$iMED %>% prcomp()
a <- hallmarkSets[hallmarkSubsets_v2] %>% unlist() %>% unique()
b <- protein_Entrez %>% filter(ENTREZID %in% a)

protein_pca2 <- proteinDat_impute$iMED %>% dplyr::select(b$OlinkID) %>% prcomp()
ggplot(data.frame(protein_pca2$x), aes(PC1, PC2, color = metadat_iMED_temp$category)) +
  geom_point(alpha=.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  ggsci::scale_color_nejm()
