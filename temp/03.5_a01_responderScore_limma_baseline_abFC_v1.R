library(limma)

get.proteinName <- function(OlinkIDs, OlinkDat) {
  # Aim: get protein names from thier Olink ID
  # input: OlinkIDs - vector of interested Olink ID protein, OlinkDat - Olink ID to protein name data 
  # output: protein names
  
  outcome <- OlinkDat %>% 
    filter(OlinkID %in% OlinkIDs)
  return(outcome$Assay)
}

# input data --------------------------
inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

# run limmar model for ab FC ===================================================
identical(metadat_healthy$name, colnames(inputDat))

designTable <- model.matrix(~ + sex + age + H1N1_abFC_combine, metadat_healthy)
res <- lmFit(inputDat, design = designTable) %>% eBayes() 

# check result
View(res$p.value)
H1N1abFC_pro <- res$p.value %>% as.data.frame() %>% filter(H1N1_abFC_combine < 0.05)
H1N1abFC_proDat <- get.proteinName(rownames(H1N1abFC_pro), cohorts_dat$proteinAnnot_all)

## Heatmap -----------------------------------
proInput <- get.datExpr(inputDat = protein_normDat,
                        metadat = metadat_healthy,
                        OlinkIDs = rownames(H1N1abFC_pro),
                        proteinAnnot = cohorts_dat$proteinAnnot_all)

heatmap_dat <- metadat_healthy %>% select(name, age_group, H1N1_reclassify) %>%
  full_join(proInput %>% t() %>% 
              as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(H1N1_reclassify, age_group) %>% 
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(H1N1_reclassify, "_", age_group)) %>% 
  column_to_rownames("name2")

heatmap_dat[, -c(1,2)] %>% t() %>% Heatmap()


## check specific proteins -----------------------------------------
proteinDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>% 
  filter(name %in% metadat$name) %>%
  column_to_rownames("name") %>% 
  select(rownames(H1N1abFC_pro))

metadat2 <- metadat_healthy %>% 
  full_join(as.data.frame(proteinDat) %>% rownames_to_column("name")) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("HL", "HH", "LL", "LH")))


proId <- "OID20462"
cowplot::plot_grid(
  ggboxplot(metadat2, x = "H1N1_reclassify", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(metadat2, x = "category", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test"),
  nrow = 2
)

# run limmar model for baseline ===================================================
identical(metadat_healthy$name, colnames(inputDat))

designTable <- model.matrix(~ + sex + age + H1N1_T1_combine, metadat_healthy)
res <- lmFit(inputDat, design = designTable) %>% eBayes() 

# check result
View(res$p.value)
H1N1_T1pro <- res$p.value %>% as.data.frame() %>% filter(H1N1_T1_combine < 0.05)
H1N1_T1proDat <- get.proteinName(rownames(H1N1_T1pro), cohorts_dat$proteinAnnot_all)

## Heatmap -----------------------------------
proInput <- get.datExpr(inputDat = protein_normDat,
                        metadat = metadat_healthy,
                        OlinkIDs = rownames(H1N1_T1pro),
                        proteinAnnot = cohorts_dat$proteinAnnot_all)

heatmap_dat <- metadat_healthy %>% select(name, age_group, H1N1_reclassify) %>%
  full_join(proInput %>% t() %>% 
              as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(H1N1_reclassify, age_group) %>% 
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(H1N1_reclassify, "_", age_group)) %>% 
  column_to_rownames("name2")

heatmap_dat[, -c(1,2)] %>% t() %>% Heatmap()

# scater plot ===========================================
avgProtein.score_abFC <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>% 
  filter(name %in% metadat$name) %>%
  column_to_rownames("name") %>% 
  select(rownames(H1N1abFC_pro)) %>%
  rowMeans(na.rm = TRUE)

avgProtein.score_T1 <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>% 
  filter(name %in% metadat$name) %>%
  column_to_rownames("name") %>% 
  select(rownames(H1N1_T1pro)) %>%
  rowMeans(na.rm = TRUE)

k <- as.data.frame(avgProtein.score_abFC) %>% rownames_to_column("name") %>%
  full_join(as.data.frame(avgProtein.score_T1) %>% rownames_to_column("name")) %>%
  full_join(metadat_healthy)

k %>% 
  ggplot(aes(x = avgProtein.score_T1, y = avgProtein.score_abFC)) + 
  geom_point(aes(col = H1N1_reclassify))+ 
  xlim(-0.5, 2) + ylim(-0.5, 2) + theme_bw()

k %>% filter(H1N1_reclassify == "LL") %>%
  ggplot(aes(x = avgProtein.score_T1, y = avgProtein.score_abFC)) + 
  geom_point(aes(col = H1N1_reclassify)) + 
  xlim(-0.5, 2) + ylim(-0.5, 2) + theme_bw()

k %>% filter(H1N1_reclassify == "LH") %>%
  ggplot(aes(x = avgProtein.score_T1, y = avgProtein.score_abFC)) + 
  geom_point(aes(col = H1N1_reclassify)) + 
  xlim(-0.5, 2) + ylim(-0.5, 2) + theme_bw()

k %>% filter(H1N1_reclassify == "HL") %>%
  ggplot(aes(x = avgProtein.score_T1, y = avgProtein.score_abFC)) + 
  geom_point(aes(col = H1N1_reclassify)) + 
  xlim(-0.5, 2) + ylim(-0.5, 2) + theme_bw()

k %>% filter(H1N1_reclassify == "HH") %>%
  ggplot(aes(x = avgProtein.score_T1, y = avgProtein.score_abFC)) + 
  geom_point(aes(col = H1N1_reclassify)) + 
  xlim(-0.5, 2) + ylim(-0.5, 2) + theme_bw()
