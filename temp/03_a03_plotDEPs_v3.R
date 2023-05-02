library(tidyverse)
library(ggpubr)
library(gplots)

# venn diagram -------------------------------------------------------------
res.venn2 <- list(
  "H1N1" = c("NBN", "SIT1", "IL15RA", "SULT2A1", "PSPN", "BTN3A2", 
             "CXADR", "CD83", "LIFR", "TNFSF12", "CXCL8", "FCRL2",
             "TNFRSF11A", "TNFRSF4", "IL12B", "CD48", "OSCAR"),
  "H3N2" = c("FABP9", "ITM2A", "FOXO1", "LSP1", "SAMD9L", "NT5C3A",
             "CLIP2", "PTX3", "HSD11B1", "LTA", "RAB6A", "KLRB1", 
             "SERPINB8", "CD79B", "TNFRSF11A", "FASLG", "CRIM1",
             "FGF19", "TPP1", "FST", "WFIKKN2"),
  "B" = c("SIT1", "TNF", "CSF3", "IL15RA", "PSPN", "BCL2L11",
          "SH2D1A", "CCL7", "LSP1", "IL7", "BTN3A2", "IL6", "FCAR",
          "CD83", "NCR1", "ICAM4", "CD4", "KRT19", "SLAMF7", "RAB6A",
          "TNFSF10", "CXCL17", "KLRB1", "SERPINB8", "CXCL8", 
          "CLEC7A", "TPSAB1", "TNFRSF11A", "VEGFA", "CCL13", "GZMA",
          "CCL11", "CXCL9", "SIGLEC1", "SPINT2", "SIRPB1", "OSCAR", "REG4", "KLRD1"))
v.table <- venn(res.venn2)


# plot ---------------------------------------------------
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name))

inputDat2 <- metadat_healthy  %>%
  left_join(inputDat)

picked_proteins <- cohorts_dat$proteinAnnot_all %>% 
  filter(Assay %in% c("NBN", "SIT1", "IL15RA", "SULT2A1", "PSPN", "BTN3A2", 
                      "CXADR", "CD83", "LIFR", "TNFSF12", "CXCL8", "FCRL2",
                      "TNFRSF11A", "TNFRSF4", "IL12B", "CD48", "OSCAR",
                      "FABP9", "ITM2A", "FOXO1", "LSP1", "SAMD9L", "NT5C3A",
                      "CLIP2", "PTX3", "HSD11B1", "LTA", "RAB6A", "KLRB1", 
                      "SERPINB8", "CD79B", "TNFRSF11A", "FASLG", "CRIM1",
                      "FGF19", "TPP1", "FST", "WFIKKN2",
                      "SIT1", "TNF", "CSF3", "IL15RA", "PSPN", "BCL2L11",
                      "SH2D1A", "CCL7", "LSP1", "IL7", "BTN3A2", "IL6", "FCAR",
                      "CD83", "NCR1", "ICAM4", "CD4", "KRT19", "SLAMF7", "RAB6A",
                      "TNFSF10", "CXCL17", "KLRB1", "SERPINB8", "CXCL8", 
                      "CLEC7A", "TPSAB1", "TNFRSF11A", "VEGFA", "CCL13", "GZMA",
                      "CCL11", "CXCL9", "SIGLEC1", "SPINT2", "SIRPB1", "OSCAR", "REG4", "KLRD1"))

plot_dat <- inputDat2 %>% 
  select(patientID, cohort, condition, age_group, category, matches("reclassify|ab"),
         picked_proteins$OlinkID)

names(plot_dat)[24:87] <- picked_proteins$Assay

# scatter plots -------------
plot_dat %>%
  ggplot(aes(x = H1N1_abFC_combine, 
             # y = TNFRSF11A
             # y = SIT1
             # y = IL15RA
             # y = OSCAR
             # y = CD83
             # y = PSPN
             # y = CXCL8
             # y = BTN3A2
             # y = LSP1
             # y = SERPINB8
             # y = TNF
             y = CD4
             )) +
  geom_point(aes(col = H1N1_reclassify),
             size = 2, alpha = 0.8,
             position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

plot_dat %>%
  ggplot(aes(x = H3N2_abFC_combine, 
             # y = TNFRSF11A
             # y = SIT1
             # y = IL15RA
             # y = OSCAR
             # y = CD83
             # y = PSPN
             # y = CXCL8
             # y = BTN3A2
             # y = LSP1
             # y = SERPINB8
             # y = TNF
             y = CD4
  )) +
  geom_point(aes(col = H1N1_reclassify),
             size = 2, alpha = 0.8,
             position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

plot_dat %>% filter(cohort == "iMED") %>%
  ggplot(aes(x = ab_B, 
             # y = TNFRSF11A
             # y = SIT1
             # y = IL15RA
             # y = OSCAR
             # y = CD83
             # y = PSPN
             # y = CXCL8
             # y = BTN3A2
             # y = LSP1
             # y = SERPINB8
             # y = TNF
             y = CD4
  )) +
  geom_point(aes(col = H1N1_reclassify),
             size = 2, alpha = 0.8,
             position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

# boxplot ---------------
my_comparision <- list(c("LL", "LH"), c("LL", "HH"), c("LH", "HH"))
strain <- "H1N1_reclassify"
strain <- "H3N2_reclassify"
strain <- "B_reclassify"

protein <- "TNFRSF11A"
protein <- "SIT1"
protein <- "IL15RA"
protein <- "OSCAR"
protein = "CD83"
protein <- "PSPN"
protein <- "CXCL8"
protein <- "BTN3A2"
protein <- "LSP1"
protein <- "SERPINB8"
protein <- "TNF"
protein <- "CD4"

plot_dat %>%
  ggboxplot(x = strain, y = protein, color = "age_group", add = "jitter") +
  stat_compare_means(comparisons = my_comparision, method = "t.test")

plot_dat %>%
  ggboxplot(x = strain, y = protein, add = "jitter") +
  stat_compare_means(comparisons = my_comparision, method = "t.test")

plot_dat %>% filter(cohort == "iMED") %>%
  ggboxplot(x = strain, y = protein, color = "age_group", add = "jitter") +
  stat_compare_means(comparisons = my_comparision, method = "t.test")


