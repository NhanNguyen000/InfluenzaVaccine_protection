rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")
#load("selected_DAMs_padj2015.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d28")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot ================================
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abT1 <- list(c("low", "high"))

mebo <- "C18H32O2"
mebo <- "C4H8O3"
mebo <- "C11H12N2O2"
mebo <- "C3H7NO2S"
mebo <- "C7H12N2O3"
mebo <- "C5H12N2O2"
mebo <- "C9H18O"
mebo <- "C11H12N2O2"
mebo <- "C3H6O4"
mebo <- "C7H10O6"
mebo <- "C9H20NO2"
mebo <- "C3H7O6P"
mebo <- "C3H7NO2S"

mebo <- selected_DAs[1]
mebo <- "C6H10O3"


mebo_classSet$`Fatty Acyls`[grepl("C20", mebo_classSet$`Fatty Acyls`)]
mebo <- "C20H38O2" # paullinic acid
mebo <- "C20H30O2" 
mebo <- "C20H36O2" 

mebo_classSet$`Fatty Acyls`[grepl("C18", mebo_classSet$`Fatty Acyls`)]
mebo <- "C18H30O2" # linolenic
mebo <- "C18H32O2" # linoleic acid
mebo <- "C18H36O2" # stearic acid

mebo_classSet$`Fatty Acyls`[grepl("C16", mebo_classSet$`Fatty Acyls`)]
mebo <- "C16H28O2"
mebo <- "C16H30O2"

mebo_classSet$`Fatty Acyls`[grepl("C14", mebo_classSet$`Fatty Acyls`)]
mebo <- "C14H26O2"
mebo <- "C14H28O2"

mebo <- "C8H14O3"

# amino acid - Arginine and proline metabolism
mebo <- "C5H9NO4"
mebo <- "C5H9NO2"
mebo <- "C5H9NO3"
mebo <- "C4H7N3O"
mebo <- "C5H7NO3"
mebo <- "C4H9NO2"

# amino acid - Valine, leucine and isoleucine degradation
mebo <- "C6H13NO2"
mebo <- "C5H8O3"
mebo <- "C4H9NO2"
mebo <- "C6H10O3"
mebo <- "C5H10O3"
mebo <- "C5H11NO2"

# amino acid - Metabolism of amino acids and derivatives
mebo <- "C10H14N5O7P"
mebo <- "C10H15N5O10P2"
mebo <- "C4H6O5"
mebo <- "C6H8N2O3"
mebo <- "C9H9NO4"
mebo <- "C6H14N2O2"
mebo <- "C5H9N3"
mebo <- "C2H7NO3S"
mebo <- "C5H7NO3"
mebo <- "C4H7NO4"
mebo <- "C2H4O3"

mebo <- "C3H7NO2S"
mebo <- "C3H6O3"

mebo <- "C5H11NO2S"

metadat_boxplot <- inputDat %>% 
  select(season, responder, c(mebo), matches("_abFC|_d0|_d28|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_d0 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

# based on abFC: NR vs R
ggboxplot(metadat_boxplot, x = "H1N1_abFC", y = mebo,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abFC, method = "t.test")

# based on abT1: low Ab vs high Ab
ggboxplot(metadat_boxplot, x = "H1N1_abBaseline", y = mebo,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abT1, method = "t.test")

# based on responders group: NR, other, TR
ggboxplot(metadat_boxplot, x = "responder", y = mebo,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_responder, method = "t.test")

# calculate the correlation --------------
# amino acid - Metabolism of amino acids and derivatives
mebo <- "C10H14N5O7P"
mebo <- "C10H15N5O10P2"
mebo <- "C4H6O5"
mebo <- "C6H8N2O3"
mebo <- "C9H9NO4"
mebo <- "C6H14N2O2"
mebo <- "C5H9N3"
mebo <- "C2H7NO3S"
mebo <- "C5H7NO3"
mebo <- "C4H7NO4"
mebo <- "C2H4O3"

mebo <- "C3H7NO2S"
cor.test(inputDat$C3H7NO2S, inputDat$C2H4O3)
cor.test(inputDat$C3H7NO2S, inputDat$C5H9N3)
cor.test(inputDat$C3H7NO2S, inputDat$C18H30O2)
cor.test(inputDat$C3H7NO2S, inputDat$C5H7NO3)
cor.test(inputDat$C3H7NO2S, inputDat$C5H11NO2S)

which(selected_DAs == "C3H7NO2S")

# calculate the ratio for FAs -------------------------------
names(mebo_Dat$ZirFlu_2019)[grep("C20", names(mebo_Dat$ZirFlu_2019))]
names(mebo_Dat$ZirFlu_2019)[grep("C18", names(mebo_Dat$ZirFlu_2019))]
names(mebo_Dat$ZirFlu_2019)[grep("C16", names(mebo_Dat$ZirFlu_2019))] # no saturated FA (C16H32O2)
names(mebo_Dat$ZirFlu_2019)[grep("C14", names(mebo_Dat$ZirFlu_2019))]


mebos <- c("C20H30O2", "C20H32O2", "C20H36O2", "C20H38O2", "C20H40O2",
           "C18H28O2", "C18H30O2", "C18H32O2", "C18H36O2",
           "C14H22O2", "C14H24O2", "C14H26O2", "C14H28O2")

metadat_boxplot <- inputDat %>% 
  select(season, responder, c(mebos), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

mebo <- "C20H30O2"
mebo <- "C20H32O2"
mebo <- "C20H36O2"
mebo <- "C20H38O2"
mebo <- "C20H40O2"

mebo <- "C18H28O2"
mebo <- "C18H30O2"
mebo <- "C18H32O2"
mebo <- "C18H36O2"

mebo <- "C14H22O2"
mebo <- "C14H24O2"
mebo <- "C14H26O2"
mebo <- "C14H28O2"
metadat_boxplot %>% ggboxplot(x = "H1N1_reclassify", y = mebo,
                              paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

a <- metadat_boxplot %>% mutate(ratio = C18H32O2 - C18H36O2)
a <- metadat_boxplot %>% mutate(ratio = C18H30O2 - C18H36O2)
a <- metadat_boxplot %>% mutate(ratio = C18H28O2 - C18H36O2)

a <- metadat_boxplot %>% mutate(ratio = C20H38O2 - C20H40O2)
a <- metadat_boxplot %>% mutate(ratio = C20H36O2 - C20H40O2)
a <- metadat_boxplot %>% mutate(ratio = C20H32O2 - C20H40O2)
a <- metadat_boxplot %>% mutate(ratio = C20H30O2 - C20H40O2)

a <- metadat_boxplot %>% mutate(ratio = C14H22O2 - C14H28O2)
a <- metadat_boxplot %>% mutate(ratio = C14H24O2 - C14H28O2)
a <- metadat_boxplot %>% mutate(ratio = C14H26O2 - C14H28O2)
a %>% ggboxplot(x = "H1N1_reclassify", y = "ratio",
                paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")


