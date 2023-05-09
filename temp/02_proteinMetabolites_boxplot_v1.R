# boxplot all_together (healthy & cirrhosis) PROTEIN DATA ------------------
# all subjects (healhty + cirrhosis)
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>% 
  column_to_rownames("name") %>% t() %>% as.data.frame %>% rownames_to_column("OlinkID") %>%
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay") %>% t() %>% as.data.frame %>% rownames_to_column("name")

metadat_boxplot <- metadat %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  left_join(inputDat)

my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("HL", "HH"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HH"))

protein <- "CD83"
ggboxplot(metadat_boxplot, x = "category", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~group) + 
 # stat_compare_means(method = "anova", label.y = 2.4) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) + # y = 3.3 or 2.8
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
# stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "B_reclassify", y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Bvictoria_reclassify", y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Byamagata_reclassify", y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

# plot protein in iMED to link with RNAseq dat
protein <- "CD83"
protein <- "IL6"

cowplot::plot_grid(
  metadat_boxplot %>% filter(cohort == "iMED") %>%
    ggboxplot(x = "H1N1_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  metadat_boxplot %>% filter(cohort == "iMED") %>%
    ggboxplot(x = "H3N2_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  metadat_boxplot %>% filter(cohort == "iMED") %>%
    ggboxplot(x = "B_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 1
)

protein <- "CD83"
protein <- "IL6"

cowplot::plot_grid(
  metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
    ggboxplot(x = "H1N1_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
    ggboxplot(x = "H3N2_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
    ggboxplot(x = "Bvictoria_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
    ggboxplot(x = "Byamagata_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 1
)

# metabolites -----------
# iMED metabolites
iMED_metabolites <- iMED$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

identical(iMED_metabolites$ionIdx, as.numeric(colnames(iMED$metabolite)))

iMED_metaboliteDat <- iMED$metabolite
colnames(iMED_metaboliteDat) <- iMED_metabolites$Formula

# ZirFlu metabolites
ZirFlu_metabolites <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

identical(ZirFlu_metabolites$ionIdx, as.numeric(colnames(ZirFlu$metabolite_dat)))

ZirFlu_metaboliteDat <- ZirFlu$metabolite_dat
colnames(ZirFlu_metaboliteDat) <- ZirFlu_metabolites$Formula

overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

metaboliteDat <- iMED_metaboliteDat %>% 
  select(overlapped_metabolites) %>% 
  as.data.frame %>% rownames_to_column("name") %>%
  full_join(ZirFlu_metaboliteDat %>% 
              select(overlapped_metabolites) %>% 
              as.data.frame %>% rownames_to_column("name"))

metadat_boxplot <- metadat %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  left_join(metaboliteDat)

metaFormula <- "C4H8O3"
metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "category", y = metaFormula ,
          paletter = "jco", add = "jitter") + facet_wrap(~group) + 
  # stat_compare_means(method = "anova", label.y = 2.4) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "H1N1_reclassify", y = metaFormula ,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) + # y = 3.3 or 2.8
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "H1N1_reclassify", y = metaFormula ,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) + # y = 3.3 or 2.8
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
