rm(list = ls())

library(tidyverse)
library(openxlsx)
# load data =======================================================================
load("/vol/projects/CIIM/Influenza/iMED/metabolic/metabolite_reclassify_timeDynamics.RData")

## sig. metabolites  
allMeboChangingSigni <- list(
  "LL_T2vsT1" = rownames(LL_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LL_T3vsT1" = rownames(LL_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LH_T2vsT1" = rownames(LH_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LH_T3vsT1" = rownames(LH_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HL_T2vsT1" = rownames(HL_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HL_T3vsT1" = rownames(HL_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HH_T2vsT1" = rownames(HH_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HH_T3vsT1" = rownames(HH_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)))
mebos <- unique(unlist(allMeboChangingSigni))

# metabolite taxonomy ================================================================
## metabolite taxonomy with HMDB ids ------------------------
mebo_taxo <- read.delim("/vol/projects/CIIM/Influenza/iMED/metabolic/db/hmdb/metabolite_HMDB_taxonomy.csv",
                        quote = "", header = TRUE)

# link the HMDB ids to Formulas
mebo_taxo_fomula <- mebo_taxo %>% rownames_to_column("CompoundID") %>% 
  inner_join(read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                       sheet = 'annotation') %>% fill(ionIdx, .direction = "down")) %>% 
  filter(Formula %in% mebos) %>%
  select(Formula, super_class, class, sub_class) %>% distinct()

length(unique(mebo_taxo_fomula$Formula)) # get all 225 formulas out of 226 formulas

mebo_taxo_fomula %>% count(super_class)
mebo_taxo_fomula %>% count(class)

## metbo stack bar charts  ------------------------------------------------
sigMebos <- allMeboChangingSigni %>%
  lapply(function(x) x %>% as.data.frame() %>% rename("valName" = ".")) %>%
  bind_rows(.id = "group")

sigMebos_taxo <- sigMebos %>% 
  left_join(mebo_taxo_fomula, by = c("valName" = "Formula"), relationship = "many-to-many")

# stack plot absolute count
sigMebos_taxo_v2 <- sigMebos_taxo %>% add_count(group, class) %>%
  mutate(metabolite_class = ifelse(n >= 10, class, "other")) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1")))

sigMebos_taxo_v2 %>%
  ggplot(aes(x = group, fill = metabolite_class)) + 
  geom_bar()

# stake plot percentage
sigMebos_taxo_v2_pct <- sigMebos_taxo_v2 %>%
  select(group, metabolite_class, valName) %>% distinct() %>%
  count(group, metabolite_class) %>%
  group_by(group) %>%
  mutate(prop = prop.table(n) * 100) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1")))

sigMebos_taxo_v2_pct %>%
  ggplot(aes(x = group, y = prop, fill = metabolite_class)) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# check the positive and negative direction separately ---------------------
allMeboChangingTstat <- list(
  "LL_T2vsT1" = LL_T2vsT1_metabol %>% select(t),
  "LL_T3vsT1" = LL_T3vsT1_metabol %>% select(t),
  "LH_T2vsT1" = LH_T2vsT1_metabol %>% select(t),
  "LH_T3vsT1" = LH_T3vsT1_metabol %>% select(t),
  "HL_T2vsT1" = HL_T2vsT1_metabol %>% select(t),
  "HL_T3vsT1" = HL_T3vsT1_metabol %>% select(t),
  "HH_T2vsT1" = HH_T2vsT1_metabol %>% select(t),
  "HH_T3vsT1" = HH_T3vsT1_metabol %>% select(t)) %>% 
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "group")

allMeboChangingTstat_selected <- allMeboChangingTstat %>% filter(valName %in% mebos) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1"))) %>%
  left_join(mebo_taxo_fomula, by = c("valName" = "Formula"), relationship = "many-to-many") %>%
  add_count(group, class) %>%
  mutate(metabolite_class = ifelse(n >= 10, class, "other"))

neg_Mebos <- allMeboChangingTstat_selected %>% filter(t < 0)
pos_Mebos <- allMeboChangingTstat_selected %>% filter(t > 0)


# stake plot percentage ---------------------------------------------
#dat_temp <- pos_Mebos
dat_temp <- neg_Mebos

plotDat_pct <- dat_temp %>%
  select(group, metabolite_class, valName) %>% distinct() %>%
  count(group, metabolite_class) %>% group_by(group) %>%
  mutate(prop = prop.table(n) * 100) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LH_T2vsT1", "HL_T2vsT1", "HH_T2vsT1",
                                   "LL_T3vsT1", "LH_T3vsT1", "HL_T3vsT1", "HH_T3vsT1")))

plotDat_pct %>% slice(grep("T2", group)) %>%
  ggplot(aes(x = group, y = prop, fill = metabolite_class)) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_pct %>% slice(grep("T3", group)) %>%
  ggplot(aes(x = group, y = prop, fill = metabolite_class)) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
