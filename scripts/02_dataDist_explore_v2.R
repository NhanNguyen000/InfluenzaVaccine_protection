# check groups
cohorts_dat$donorInfo_all %>% count(cohort, season, disease, age_group)

# NEW CODE ---------------------------
library(ggpubr)

metadat <- cohorts_dat$donorInfo_all %>%
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# age distribution ----------------------------
plot_dat <- cohorts_dat$donorInfo_all %>%
  mutate(group = paste(cohort, season, disease, sep = "_"))
ggplot(plot_dat, aes(x=age, fill=group)) + geom_density(alpha=.3) + xlim(20, 81)
ggplot(plot_dat %>% filter(season %in% c(NA, "2019")), 
       aes(x=age, fill=group)) + geom_density(alpha=.3)

ggplot(plot_dat, aes(x=age, fill=cohort)) + geom_density(alpha=.3)+ xlim(20, 81)

plot_dat %>% filter(disease == "healthy") %>% 
  ggplot(aes(x=age, fill=cohort)) + geom_density(alpha=.3) + xlim(20, 81)

plot_dat %>% 
  mutate(cohort2 = ifelse(cohort == "iMED", cohort, paste0(cohort, "_", disease))) %>% 
  ggplot(aes(x=age, fill=cohort2)) + geom_density(alpha=.3) + xlim(20, 81)

plot_dat %>% 
  mutate(cohort2 = ifelse(cohort == "iMED", cohort, paste0(cohort, "_", disease))) %>% 
  filter(cohort2 != "ZirFlu_cirrhosis") %>%
  ggplot(aes(x=age, fill=cohort2)) + geom_density(alpha=.3) + xlim(20, 81)

## sex and age ----------------------------

metadat %>% 
  ggboxplot(x = "sex", y = "age",
          paletter = "jco", add = "jitter") + 
  facet_wrap(~group) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

metadat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "sex", y = "age",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~group) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

metadat %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "sex", y = "age",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~group) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

## baseline ----------------------------
strain <- "H1N1_T1_combine"
strain <- "H3N2_T1_combine" 

# with age
metadat %>% 
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("old", "young"), c("young", "cirrhosis"), c("old", "cirrhosis")),
                     method = "t.test")

metadat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7)))

metadat %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("old", "young"), c("young", "cirrhosis"), c("old", "cirrhosis")),
                     method = "t.test")
# with age + sex
metadat %>% 
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

metadat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

metadat %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(ab) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

## abFC ------------------------------
strain <- "H1N1_abFC_combine"
strain <- "H3N2_abFC_combine"

# with age
metadat %>% 
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("old", "young"), c("young", "cirrhosis"), c("old", "cirrhosis")),
                     method = "t.test")

metadat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7)))

metadat %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("old", "young"), c("young", "cirrhosis"), c("old", "cirrhosis")),
                     method = "t.test")

metadat %>% 
  ggboxplot(x = "group", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~cohort) +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("old", "young"), c("young", "cirrhosis"), c("old", "cirrhosis")),
                     method = "t.test")

# with age + sex
metadat %>% 
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

metadat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

metadat %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "sex", y = strain,
            paletter = "jco", add = "jitter") +
  facet_wrap(~group) +
  ylab(paste0("log2(abFC) - ", substr(strain, 1, 7))) + 
  stat_compare_means(comparisons = list(c("f", "m")),
                     method = "t.test")

## abFC vs. ab baselines ----------
metadat %>%
  ggplot(aes(x = H1N1_abFC_combine, y = H1N1_T1_combine)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  xlab("log2(abFC) - H1N1") + ylab("log2(ab) - H1N1") + 
  theme_classic()

# H1N1 per cohort
metadat %>% filter(cohort == "iMED") %>%
  ggplot(aes(x = H1N1_abFC_combine, y = H1N1_T1_combine)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  xlab("log2(abFC) - H1N1") + ylab("log2(ab baseline) - H1N1") + 
  theme_classic()

metadat %>% filter(cohort == "ZirFlu") %>%
  ggplot(aes(x = H1N1_abFC_combine, y = H1N1_T1_combine)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  xlab("log2(abFC) - H1N1") + ylab("log2(ab baseline) - H1N1") + 
  theme_classic()

# B per cohort
metadat %>% filter(cohort == "iMED") %>%
  ggplot(aes(x = iMED_B_T1, y = ab_B)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  xlab("log2(abFC) - B") + ylab("log2(ab baseline) - B") + 
  theme_classic()

metadat %>% filter(cohort == "ZirFlu") %>%
  ggplot(aes(x = Bvictoria_abFC, y = Bvictoria_T1)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  xlab("log2(abFC) - Bvictoria") + ylab("log2(ab baseline) - Bvictoria") + 
  theme_classic()

# classification (responder vs. non-responder) ------------------
# iMED - H1N1
plotDat <- metadat %>% filter(cohort == "iMED") %>%
  dplyr::select(patientID, cohort, group, iMED_H1N1_T1, iMED_H1N1_T4, ab_H1N1) %>%
  mutate(category = ifelse(ab_H1N1 >= 4, "response", "non-response")) %>%
  select(-ab_H1N1) %>%
  gather(matches("_T"), key = time, value = HAItiter) %>%
  arrange(patientID) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("iMED_H1N1_", "", time))

# ZirFlu - H1N1
plotDat <- metadat %>% filter(cohort == "ZirFlu") %>%
  dplyr::select(patientID, cohort, group, H1N1_T1, H1N1_T2, H1N1_abFC) %>%
  mutate(category = ifelse(H1N1_abFC >= 4, "response", "non-response")) %>%
  select(-H1N1_abFC) %>%
  gather(matches("_T"), key = time, value = HAItiter) %>%
  arrange(patientID) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("H1N1_", "", time))

# plot dat
plotDat %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(category), nrow = 1) +
  ylab("HAI titer - H1N1")

# reclassification --------------------------
# iMED - H1N1
plotDat <- metadat %>% filter(cohort == "iMED") %>%
  dplyr::select(patientID, cohort, group, iMED_H1N1_T1, iMED_H1N1_T4, H1N1_reclassify) %>%
  gather(matches("_T"), key = time, value = HAItiter) %>%
  arrange(patientID) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("iMED_H1N1_", "", time))

# ZirFlu - H1N1
plotDat <- metadat %>% filter(cohort == "ZirFlu") %>%
  dplyr::select(patientID, cohort, group, H1N1_T1, H1N1_T2, H1N1_reclassify) %>%
  gather(matches("_T"), key = time, value = HAItiter) %>%
  arrange(patientID) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("H1N1_", "", time))

# plot data
plotDat %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(reclassify), nrow = 1) +
  ylab("HAI titer - H1N1")

# PCA for protein & metabolites ---------------------
library(caret)
library(tidyverse)

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
         color = groupType) +
    stat_ellipse()
  
}
## make PCA for protein ------------------
# impute protein
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>% 
  arrange(factor(name, levels = metadat$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay") %>% drop_na()

get.pca_plot(pca_dat, metadata, "condition")

iMED_protein <- list() 
iMED_protein$metadata <- metadat %>% filter(cohort == "iMED")
iMED_protein$dat <- inputDat %>% select(iMED_protein$metadata$name)

a <- prcomp(iMED_protein$dat)

pca_dat <- predict(preProcess(x=iMED_protein$dat, method = "knnImpute"),
                       iMED_protein$dat) %>% prcomp()
summary(pca_dat)$importance[1:3, 1:5]
get.pca_plot(a, iMED_protein$metadata, "category")

metabolite_pca <- ZirFlu$metabolite_dat %>% prcomp()
summary(metabolite_pca)$importance[1:3, 1:5]
get.cowplot_pca(pca_dat = metabolite_pca, metadata)
