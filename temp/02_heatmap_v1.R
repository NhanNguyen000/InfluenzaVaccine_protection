# load libraries & functions --------------------------
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
# annotation color -------------------------------------------------------------
col_cohort <- c("iMED" = "blue", "ZirFlu" = "red")
col_age <- c("young" = "blue", "old" = "red")
col_condition <- c("healthy" = "orange", 
                   "compensated cirrhosis" = "darkblue", 
                   "decompensated cirrhosis"  = "darkred")
col_reclasify <- c("LL" = "black", "LH" = "green", "HL" = "blue", "HH" = "red")
col_category <- c("NR" = "grey90", "Other" = "khaki3", "TR" = "peru")
col_inflam <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

# input data (no NA) & metadata -------------------------------------------------
## ZirFlu 2019 -----------------------------------
inputDat <- proteinDat_impute$ZirFlu_2019 %>% 
  t() %>% as_tibble()

metadat <- cohorts_dat$donorSample_all %>% filter(time == "T1") %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(name %in% colnames(inputDat)) %>% drop_na(category) %>%
  left_join(inflamScore_dat)

inputDat <- inputDat %>% as_tibble() %>% dplyr::select(metadat$name)

## all ZirFlu & iMED -----------------------------
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  column_to_rownames("name") %>% t() %>% 
  as_tibble() %>% drop_na()

inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  column_to_rownames("name") %>% 
  dplyr::select(inflam.Proteins$OlinkID) %>% 
  t() %>% as_tibble() %>% drop_na()

metadat <- cohorts_dat$donorSample_all %>% filter(time == "T1") %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(name %in% colnames(inputDat)) %>% drop_na(category) %>%
  left_join(inflamScore_dat)

inputDat <- inputDat %>% as_tibble() %>% dplyr::select(metadat$name)

## all ZirFlu & iMED - only healthy -----------------------------
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  column_to_rownames("name") %>% 
  t() %>% as_tibble(rownames = NA) %>% 
  rownames_to_column("name") %>% drop_na() %>% column_to_rownames("name")

inflamScore_dat <- inflamScore %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) 

metadat_healthy <- cohorts_dat$donorSample_all %>% filter(time == "T1") %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(name %in% colnames(inputDat)) %>% drop_na(category) %>%
  filter(disease == "healthy") %>%
  left_join(inflamScore_dat) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate(age_group2 = ifelse(cohort == "iMED", "iMED_>64",
                             ifelse(age >= 64, "ZirFlu_>64",
                                    ifelse(age >= 40, "ZirFlu_40-64", "ZirFlu_20-40")))) %>%
  mutate(age_group2 = factor(age_group2, levels = c("ZirFlu_20-40", "ZirFlu_40-64", "ZirFlu_>64", "iMED_>64")))

inputDat <- inputDat %>% dplyr::select(metadat$name)

## all ZirFlu & iMED - only healthy subjects & inflammation proteins -----------------------------
load("gsea_hallmarktSet.RData") # load hallmark gene sets (prepare from from MSigDB) & convert protein - Entrez gene ID
inflam.Proteins <- protein_Entrez %>% 
  filter(ENTREZID %in% c(hallmarkSets[hallmarkSubsets_v2] %>% as_vector() %>% unique()))

inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  column_to_rownames("name") %>%
  dplyr::select(inflam.Proteins$OlinkID) %>% 
  t() %>% as_tibble() %>% drop_na()

metadat <- cohorts_dat$donorSample_all %>% filter(time == "T1") %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(name %in% colnames(inputDat)) %>% drop_na(category) %>%
  filter(disease == "healthy")

inputDat <- inputDat %>% as_tibble() %>% dplyr::select(metadat_healthy$name)

inflam_aveProteins <- inputDat %>% colMeans()
inflam_aveProteinDat <- metadat %>% 
  full_join(as.data.frame(inflam_aveProteins) %>% rownames_to_column("name"))

metadat_healthy <- metadat_healthy %>% 
  left_join(inflam_aveProteinDat %>% select(name, inflam_aveProteins))

## all ZirFlu & iMED - only healthy subjects & inflammation leading proteins -----------------------------
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  column_to_rownames("name") %>%
  dplyr::select(leadingProteins_TRvsNR_dat$OlinkID) %>% 
  t() %>% as_tibble() %>% drop_na()

inflam_aveLeadingProteins <- inputDat %>% colMeans()
inflam_aveLeadingProteinDat <- metadat_healthy %>% 
  full_join(as.data.frame(inflam_aveLeadingProteins) %>% rownames_to_column("name"))

metadat_healthy <- metadat_healthy %>% 
  left_join(inflam_aveLeadingProteinDat %>% select(name, inflam_aveLeadingProteins))

# correct for the variants
aveProteinExpr <- inputDat %>% rowMeans()
inputDat_subtracAvg <- inputDat - aveProteinExpr
inflam_aveLeadingProteins_subtracAvg <- inputDat_subtracAvg %>% colMeans()

inflam_aveLeadingProteins_subtracAvgDat <- metadat_healthy %>% 
  full_join(as.data.frame(inflam_aveLeadingProteins_subtracAvg) %>% rownames_to_column("name"))

metadat_healthy <- metadat_healthy %>% 
  left_join(inflam_aveLeadingProteins_subtracAvgDat %>% select(name, inflam_aveLeadingProteins_subtracAvg))

# heatmap ----------------------------------------------------------------------
identical(metadat$name, colnames(inputDat)) # check the sampleorder = TRUE

column_ha = HeatmapAnnotation(
  # cohort = metadat$cohort,
  age_group = metadat$age_group,
  condition = metadat$condition,
  # H1N1 = metadat$H1N1$reclassify,
  # H3N2 = metadat$H3N2$reclassify,
  # Bvic = metadat$Bvictoria$reclassify,
  # Byam = metadat$Byamagata$reclassify,
  category = metadat$category,
  inflamScore = metadat$inflamScore.hallmarkSubsets_v2,
  #avg_inflam.Protein = inflam_aveProteins,
  col = list(
    # cohort = col_cohort,
    age_group = col_age,
    condition = col_condition,
    # H1N1 = col_reclasify, 
    # H3N2 = col_reclasify, 
    # Bvic = col_reclasify, 
    # Byam = col_reclasify,
    category = col_category,
    inflamScore = col_inflam,
    inflamProtein = col_inflam))
Heatmap(inputDat,
        heatmap_legend_param = list(at = c(-5, 0, 5)))

Heatmap(inputDat,
        heatmap_legend_param = list(at = c(-5, 0, 5)),
       # show_row_names = FALSE, 
        top_annotation = column_ha )

Heatmap(inputDat,
        heatmap_legend_param = list(at = c(-5, 0, 5)),
        top_annotation = column_ha )

# check the distribution
ggplot(metadat, aes(x=inflamScore.hallmarkSubsets_v2)) + 
  geom_density(alpha=.3) +
  geom_bar()

ggplot(metadat, aes(x = age, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_point(aes(col = category))

metadat %>% filter(category != "Other") %>% 
  ggplot(aes(x = age, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_point(aes(col = category)) + theme_bw()

ggplot(inflam_aveProteinDat, aes(x=inflam_aveProteins)) + 
  geom_density(alpha=.3) +
  geom_bar()

inflam_aveProteinDat %>% filter(category != "Other") %>% 
  ggplot(aes(x = age, y = inflam_aveProteins)) + 
  geom_point(aes(col = category)) + theme_bw()

# split inflamScore into 3 group
medata2 <- metadat %>%
  mutate(inflam_group1 = as.numeric(cut_number(metadat$inflamScore.hallmarkSubsets_v2, 3))) %>%
  mutate(inflam_group2 = ifelse(inflam_group1 == 1, "inflam.low",
                                ifelse(inflam_group1 ==2, "inflam.mid", "inflam.high"))) %>%
  mutate(inflam_group2 = factor(inflam_group2, levels = c("inflam.low", "inflam.mid", "inflam.high")))

# check abFC
medata2 %>% ggplot(aes(x = inflam_group2, y = H1N1_abFC_combine)) + geom_boxplot()
boxplot(medata2$H1N1_abFC_combine ~ medata2$inflam_group2)
boxplot(medata2$H1N1_abFC_combine ~ medata2$inflam_group2 + medata2$age_group)

medata2 %>% filter(age_group == "young") %>%
  ggplot(aes(x = inflam_group2, y = H1N1_abFC_combine)) + 
  geom_boxplot() + theme_bw()

medata2 %>% 
  ggplot(aes(x = inflamScore.hallmarkSubsets_v2, y = H1N1_abFC_combine)) + 
  geom_point(aes(col = category)) + theme_bw()
medata2 %>% 
  ggplot(aes(x = inflamScore.hallmarkSubsets_v2, y = H3N2_abFC_combine)) + 
  geom_point(aes(col = category)) + theme_bw()

# all samples ------------------
boxplot(metadat$inflamScore.hallmarkSubsets_v2 ~ metadat$age)
boxplot(metadat$inflamScore.hallmarkSubsets_v2 ~ metadat$age_group)
boxplot(metadat$inflamScore.hallmarkSubsets_v2 ~ metadat$disease + metadat$age_group)
boxplot(metadat$inflamScore.hallmarkSubsets_v2 ~ metadat$category + metadat$disease)

metadat2 <- metadat %>% filter(age_group == "young")
boxplot(metadat2$inflamScore.hallmarkSubsets_v2 ~ metadat2$category + metadat2$disease)

boxplot(metadat$inflamScore.hallmarkSubsets_v2 ~ metadat$category + metadat$disease)

metadat %>% 
  ggplot(aes(x = disease, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_boxplot(aes(col = category)) + theme_bw() + facet_wrap(~age_group)

Heatmap(inputDat %>% t())

plot_dat <- metadat %>% mutate(group2 = paste0(age_group, "_", disease)) %>%
  left_join(proteinDat_impute %>% 
              lapply(function(x) x %>% rownames_to_column("name")) %>% 
              purrr::reduce(full_join))

plot_dat %>%
  ggplot(aes(x = H1N1_abFC_combine, y = OID20611)) +
  geom_point()

library(ggpubr)
plot_dat %>% filter(age_group == "young") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") # TNFSF10
plot_dat %>% 
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", add = "reg.line") +
  stat_cor(method = "pearson")

plot_dat %>% filter(category == "TR") %>%
  ggscatter(x = "H3N2_abFC_combine", y = "OID20611", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") 

plot_dat %>% 
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson")  # TNFSF10

plot_dat %>% filter(disease == "healthy") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson")  # TNFSF10

plot_dat %>% filter(age_group == "young") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") # TNFSF10

plot_dat %>% filter(category == "TR") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") 

plot_dat %>% filter(age_group == "young") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "category", add = "reg.line") +
  stat_cor(aes(color = category), method = "pearson") 

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20611", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson") 


plot_dat %>% 
  ggscatter(x = "H1N1_abFC_combine", y = "OID20563", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") # IL6

plot_dat %>% filter(disease == "healthy") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20563", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") # IL6

plot_dat %>% filter(age_group == "young") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20563", 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") # IL6

plot_dat %>% filter(category == "TR" & disease == "healthy") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20563", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson")

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20563", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson") 

plot_dat %>% 
  ggscatter(x = "H1N1_abFC_combine", y = "OID20562", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") # IL15

plot_dat %>% filter(disease == "healthy") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20562", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") 

plot_dat %>% filter(age_group == "young") %>%
  ggscatter(x = "H3N2_abFC_combine", y = "OID20562", 
            color = "age_group", add = "reg.line") +
  stat_cor(aes(color = age_group), method = "pearson") 

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20562", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson") 

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20624", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson")  # TNFSF12

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20631", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson") # CXCL8

plot_dat %>% filter(group2 != "old_cirrhosis") %>%
  ggscatter(x = "H1N1_abFC_combine", y = "OID20477", 
            color = "group2", add = "reg.line") +
  stat_cor(aes(color = group2), method = "pearson") # IL17C


# BOXPLOT & test ---------------
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )

# young
metadat2 <- metadat_healthy %>% 
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  filter(age_group == "young" & disease == "healthy")

ggboxplot(metadat2, x = "category", y = "inflamScore.hallmarkSubsets_v2",
          color = "category", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(metadat2, x = "category", y = "inflamScore.hallmarkSubsets_v2",
          color = "category", palette = "jco") +
  stat_compare_means(method = "anova")

# old
metadat2 <- metadat %>% 
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  filter(age_group == "old" & disease == "healthy")

ggboxplot(metadat2, x = "category", y = "inflamScore.hallmarkSubsets_v2",
          color = "category", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(metadat2, x = "category", y = "inflamScore.hallmarkSubsets_v2",
          color = "category", palette = "jco") +
  stat_compare_means(method = "anova")
# combine young & old

ggboxplot(metadat_healthy, 
          x = "category", 
          y = "inflamScore.hallmarkSubsets_v2",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(metadat_healthy, 
          x = "category", 
          y = "inflamScore.hallmarkSubsets_v2",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_healthy, 
          x = "age_group2", 
          y = "inflamScore.hallmarkSubsets_v2",
          color = "category",
          paletter = "jco", add = "jitter") +
  stat_compare_means(aes(group = category))

ggboxplot(metadat_healthy, 
          x = "category", 
          y = "inflam_aveProteins",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_healthy, 
          x = "age_group2", 
          y = "inflam_aveProteins",
          color = "category",
          paletter = "jco", add = "jitter") +
  stat_compare_means(aes(group = category))

ggboxplot(metadat_healthy, 
          x = "category", 
          y = "inflam_aveLeadingProteins",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_healthy, 
          x = "age_group2", 
          y = "inflam_aveLeadingProteins",
          color = "category",
          paletter = "jco", add = "jitter") +
  stat_compare_means(aes(group = category))

ggboxplot(metadat_healthy, 
          x = "category", 
          y = "inflam_aveLeadingProteins_subtracAvg",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_healthy, 
          x = "age_group2", 
          y = "inflam_aveLeadingProteins_subtracAvg",
          color = "category",
          paletter = "jco", add = "jitter") +
  stat_compare_means(aes(group = category))

# only TR
metadat2 <- metadat %>% filter(category == "TR") %>% 
  mutate(group2 = paste0(age_group, "_", disease))

ggboxplot(metadat2, x = "group2", y = "inflamScore.hallmarkSubsets_v2",
          color = "group2", palette = "jco") +
  stat_compare_means(method = "anova")

my_comparisons <- list( c("young_healthy", "old_healthy"), 
                        c("young_healthy", "young_cirrhosis"),
                        c("old_healthy", "old_cirrhosis"))
ggboxplot(metadat2, x = "group2", y = "inflamScore.hallmarkSubsets_v2",
          color = "group2", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons)

# heatmap ----------------------------------------------------------------------
metadat_inflamScore <- list()
metadat_inflamScore$ZirFlu_2019 <- add.inflamScore(metadat = metadat_sampleT1$ZirFlu_2019, 
                                                   sample_col = "probenID", 
                                                   inflamScore = inflamScore$ZirFlu_2019)
metadat_inflamScore$ZirFlu_2020 <- add.inflamScore(metadat = metadat_sampleT1$ZirFlu_2020, 
                                                   sample_col = "probenID", 
                                                   inflamScore = inflamScore$ZirFlu_2020)
metadat_inflamScore$iMED <- add.inflamScore(metadat = metadat_sampleT1$iMED, 
                                            sample_col = "name", 
                                            inflamScore = inflamScore$iMED)

# heatmap - grop color
col_condition <- c("healthy" = "orange", 
                   "compensated cirrhosis" = "darkblue", 
                   "decompensated cirrhosis"  = "darkred")
col_reclasify <- c("LL" = "black", "LH" = "green", "HL" = "blue", "HH" = "red")
col_category <- c("NR" = "grey90", "Other" = "khaki3", "other" = "khaki3", "TR" = "peru")
col_inflam <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

# ZirFlu 2019
column_ha = HeatmapAnnotation(
  condition = metadat_inflamScore$ZirFlu_2019$condition,
  H1N1 = metadat_inflamScore$ZirFlu_2019$H1N1_reclassify,
  H3N2 = metadat_inflamScore$ZirFlu_2019$H3N2_reclassify,
  Bvic = metadat_inflamScore$ZirFlu_2019$Bvictoria_reclassify,
  Byam = metadat_inflamScore$ZirFlu_2019$Byamagata_reclassify,
  category = metadat_inflamScore$ZirFlu_2019$category,
  inflam.score = metadat_inflamScore$ZirFlu_2019$inflamScore.hallmarkSubsets_v2,
  col = list(condition = col_condition,
             H1N1 = col_reclasify, 
             H3N2 = col_reclasify, 
             Bvic = col_reclasify, 
             Byam = col_reclasify,
             category = col_category,
             inflam.score = col_inflam))
Heatmap(fgseaResTab$ZirFlu_2019$hallmarkSets,
        top_annotation = column_ha )

# ZirFlu 2020
column_ha = HeatmapAnnotation(
  condition = metadat_inflamScore$ZirFlu_2020$condition,
  H1N1 = metadat_inflamScore$ZirFlu_2020$H1N1_reclassify,
  H3N2 = metadat_inflamScore$ZirFlu_2020$H3N2_reclassify,
  Bvic = metadat_inflamScore$ZirFlu_2020$Bvictoria_reclassify,
  Byam = metadat_inflamScore$ZirFlu_2020$Byamagata_reclassify,
  category = metadat_inflamScore$ZirFlu_2020$category,
  inflam.score = metadat_inflamScore$ZirFlu_2020$inflamScore.hallmarkSubsets_v2,
  col = list(condition = col_condition,
             H1N1 = col_reclasify, 
             H3N2 = col_reclasify, 
             Bvic = col_reclasify, 
             Byam = col_reclasify,
             category = col_category,
             inflam.score = col_inflam))
Heatmap(fgseaResTab$ZirFlu_2020$hallmarkSets,
        top_annotation = column_ha )

# iMED
column_ha = HeatmapAnnotation(
  condition = metadat_inflamScore$iMED$condition,
  H1N1 = metadat_inflamScore$iMED$H1N1_reclassify,
  H3N2 = metadat_inflamScore$iMED$H3N2_reclassify,
  B = metadat_inflamScore$iMED$B_reclassify,
  category = metadat_inflamScore$iMED$responder,
  inflam.score = metadat_inflamScore$iMED$inflamScore.hallmarkSubsets_v2,
  col = list(condition = col_condition,
             H1N1 = col_reclasify, 
             H3N2 = col_reclasify, 
             B = col_reclasify, 
             category = col_category,
             inflam.score = col_inflam))
Heatmap(fgseaResTab$iMED$hallmarkSets,
        top_annotation = column_ha )

