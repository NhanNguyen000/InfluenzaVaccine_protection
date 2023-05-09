# load libraries & functions --------------------------
library(ggpubr)

get.ggboxplot <- function(dat, x_col, y_col) {
  # Description: make a boxplot
  
  # Arguments: 
  # dat - data with patients in row, proteins in column
  # x_col - name of the colum with group/class information (x axis)
  # y_col - name of the colum with numeric value - will show in boxplot (y axis)
  
  ggboxplot(dat, x = x_col, y = y_col,
            paletter = "jco", add = "jitter")
}

add.inflamScore <- function(metadat, sample_col, inflamScore) {
  # Description: add the inflamation score to the metadata
  
  # Arguments: 
  # metadat - table with patient information
  # sample_col - name of the column which matched with the sample id in the inflamScore data
  # inflamScore - table with inflamation socre in column, and sample id in rownames 
  
  # Returns: outcome - metadata with added inflamation score
  
  outcome <- metadat %>% mutate_at(sample_col, as.character) %>%
    full_join(inflamScore %>% rownames_to_column(sample_col))
  return(outcome)
}
# prepare data ----------------------------------------------------------------------------
load("inflamScore_NESfgsae.RData")

inflamScore_allCohort <- cohorts_dat$donorInfo_all %>% 
  left_join(inflamScore.NESfgsea) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("HL", "HH", "LL", "LH")))

# Density plot for the inflammation score --------------------------------------------
inflamScore_allCohort %>%
  ggplot(aes(x=inflamScore.hallmarkSubsets_v2, fill=cohort)) + geom_density(alpha=.3)

inflamScore_allCohort %>% filter(cohort == "ZirFlu") %>%
  ggplot(aes(x=inflamScore.hallmarkSubsets_v2, fill= season)) + geom_density(alpha=.3)

# boxplot across age 20-81 - all subjects
boxplot(inflamScore.hallmarkSubsets_v2 ~ age,
        data = inflamScore_allCohort)

inflamScore_allCohort %>%
  ggplot(aes(x = disease, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_boxplot(aes(col = age_group)) + theme_bw()

# plot for 3 categories: NR, Other, TR --------------------------------------------
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )

# boxplot across age 20-81 - all healthy subjects
inflamScore_allCohort_healthy <- inflamScore_allCohort %>% filter(disease == "healthy")

boxplot(inflamScore.hallmarkSubsets_v2 ~ age,
        data = inflamScore_allCohort_healthy)

inflamScore_allCohort_healthy %>%
  ggplot(aes(x = category, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_boxplot() + theme_bw()

inflamScore_allCohort_healthy %>%
  ggplot(aes(x = age_group, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_boxplot(aes(col = category)) + theme_bw()

# boxplot across age 20-81 - all cirrhosis subjects
inflamScore_allCohort_cirrhosis <- inflamScore_allCohort %>% filter(disease == "cirrhosis")

boxplot(inflamScore.hallmarkSubsets_v2 ~age,
        data = inflamScore_allCohort_cirrhosis)

inflamScore_allCohort_cirrhosis %>%
  ggplot(aes(x = age_group, y = inflamScore.hallmarkSubsets_v2)) +
  geom_boxplot() + theme_bw()

inflamScore_allCohort_cirrhosis %>%
  ggplot(aes(x = age_group, y = inflamScore.hallmarkSubsets_v2)) +
  geom_boxplot(aes(col = category)) + theme_bw()

# with statistic test
inflamScore_allCohort_healthy <- inflamScore_allCohort %>% filter(disease == "healthy")


ggboxplot(inflamScore_allCohort_healthy %>% filter(disease == "healthy"), 
          x = "category", 
          y = "inflamScore.hallmarkSubsets_v2",
          color = "category",
          paletter = "jco", add = "jitter",
          facet.by = "age_group") +
  stat_compare_means(comparisions = my_comparisons, method = "t.test")

ggplot(inflamScore_allCohort_healthy, aes(x = category, y = inflamScore.hallmarkSubsets_v2)) +
  facet_wrap(~age_group) + 
  stat_compare_means(comparisions = my_comparisons)


# inflam_avg protein expression ------------------------
boxplot(inflam_aveProteinDat$inflam_aveProteins ~ inflam_aveProteinDat$age)
boxplot(inflam_aveProteinDat$inflam_aveProteins ~ inflam_aveProteinDat$age_group)
boxplot(inflam_aveProteinDat$inflam_aveProteins ~ inflam_aveProteinDat$category + inflam_aveProteinDat$age_group)

# scater plot plot with age ------------------------------------
inflamScore_allCohort %>% 
  ggplot(aes(x = age, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_point(aes(col = category))

inflamScore_allCohort %>% filter(disease == "healthy") %>%
  ggplot(aes(x = age, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_point(aes(col = category))

inflamScore_allCohort %>% filter(disease != "healthy") %>%
  ggplot(aes(x = age, y = inflamScore.hallmarkSubsets_v2)) + 
  geom_point(aes(col = category))

inflamScore_proteinDat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$HAItiter_all) %>%
  left_join(proteinDat %>% 
              lapply(function(x) x %>% rownames_to_column("name")) %>% 
              purrr::reduce(full_join) %>%
              left_join(cohorts_dat$donorSample_all))

inflamScore_proteinDat_v2 <- inflamScore_proteinDat %>% filter(disease == "healthy")
inflamScore_proteinDat_v2 %>% 
  ggplot(aes(x = age, y = OID20562)) + 
  geom_point(aes(col = category))
boxplot(OID20562 ~ age, data = inflamScore_proteinDat_v2)
boxplot(OID20562 ~ category + age_group, data = inflamScore_proteinDat_v2) # IL15

inflamScore_proteinDat_v2 %>% 
  ggplot(aes(x = age, y = OID20611)) + 
  geom_point(aes(col = category)) # TNFSF10
boxplot(OID20611 ~ age, data = inflamScore_proteinDat_v2)
boxplot(OID20611 ~ category + age_group, data = inflamScore_proteinDat_v2) 

inflamScore_proteinDat_v2 %>% 
  ggplot(aes(x = age, y = OID20624)) + 
  geom_point(aes(col = category)) # TNFSF12
boxplot(OID20624 ~ age, data = inflamScore_proteinDat_v2)
boxplot(OID20624 ~ category + age_group, data = inflamScore_proteinDat_v2) 


inflamScore_proteinDat_v2 %>% 
  ggplot(aes(x = age, y = OID20631)) + 
  geom_point(aes(col = category)) # CXCL8
boxplot(OID20631 ~ age, data = inflamScore_proteinDat_v2)
boxplot(OID20631 ~ category + age_group, data = inflamScore_proteinDat_v2) 


# boxplot for inflamation score for 4 reclassification ----------------------------------------
inflamScore_temp <- "inflamScore.hallmarkSubsets_v2"
my_comparisons_v2 <- list( c("HH", "HL"), c("LL", "LH") )

plotDat <- inflamScore_allCohort
plotDat <- inflamScore_allCohort_healthy
cowplot::plot_grid(
  get.ggboxplot(plotDat, x_col = "H1N1_reclassify", y_col = inflamScore_temp),
  get.ggboxplot(plotDat, x_col = "H3N2_reclassify", y_col = inflamScore_temp),
  get.ggboxplot(plotDat, x_col = "Bvictoria_reclassify2", y_col = inflamScore_temp),
  get.ggboxplot(plotDat, x_col = "Byamagata_reclassify2", y_col = inflamScore_temp),
  nrow = 1
)

cowplot::plot_grid(
  ggboxplot(plotDat, x = "H1N1_reclassify", y = inflamScore_temp,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(plotDat, x = "H3N2_reclassify", y = inflamScore_temp,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group )+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(plotDat, x = "Bvictoria_reclassify2", y = inflamScore_temp,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(plotDat, x = "Byamagata_reclassify2", y = inflamScore_temp,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 2
)

plotDat_ZirFlu <- inflamScore_allCohort_healthy %>% filter(cohort == "ZirFlu")
cowplot::plot_grid(
  get.ggboxplot(plotDat_ZirFlu, x_col = "H1N1_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  get.ggboxplot(plotDat_ZirFlu, x_col = "H3N2_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  get.ggboxplot(plotDat_ZirFlu, x_col = "Bvictoria_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  get.ggboxplot(plotDat_ZirFlu, x_col = "Byamagata_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 1
)

plotDat_iMED <- inflamScore_allCohort_healthy %>% filter(cohort == "iMED")
cowplot::plot_grid(
  get.ggboxplot(plotDat_iMED, x_col = "H1N1_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  get.ggboxplot(plotDat_iMED, x_col = "H3N2_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  get.ggboxplot(plotDat_iMED, x_col = "B_reclassify", y_col = inflamScore_temp)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 1
)

ggboxplot(plotDat, x = "H1N1_reclassify", y = inflamScore_temp,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2)

# check the cor(inflamScore, HAI reponse) ---------------------------------------
cor.test(inflamScore_allCohort$inflamScore.hallmarkSubsets_v2, 
         inflamScore_allCohort$H1N1_abFC_combine)

