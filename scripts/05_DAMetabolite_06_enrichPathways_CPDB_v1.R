rm(list = ls())
library(tidyverse)

# load the data ----------------
ora_pathways <- read.table("processedDat/ORA_results_DAM_4reclass_to2groups.tab", sep = "\t", header = TRUE) %>% filter(q.value < 0.05) 
ora_pathways <- read.table("processedDat/ORA_results_DAM_4groups.tab", sep = "\t", header = TRUE) %>% filter(q.value < 0.05) 

# plot the sig. pathway (padj < 0.05) -----------------
ora_pathways %>% 
  ggplot(aes(x = -log10(q.value), y = pathway, size = `effective_size`)) + 
  geom_point() + xlim(c(0, 3.5)) + theme_classic()

# check the pathway ---------------------
# option 1 group pathways into bigser group
pathways <- list()
pathways$fatty_acids <- ora_pathways %>% slice(grep("fatty acid", pathway))
pathways$linolenic <- ora_pathways %>% slice(grep("Linolenic", pathway))
pathways$phospholipases <- ora_pathways %>% slice(grep("phospholipases", pathway))
pathways$amino_acids <- ora_pathways %>%  slice(grep("amino acids", pathway))
pathways$valine <- ora_pathways %>% slice(grep("Valine", pathway))
pathways$arginine <- ora_pathways %>% slice(grep("Arginine", pathway))
pathways$sphingo <- ora_pathways %>% slice(grep("sphingo", pathway))
pathways$gAlpha <- ora_pathways %>% slice(grep("G alpha", pathway))
pathways$retinol <- ora_pathways %>% slice(grep("retinol", pathway))
pathways$steroid <- ora_pathways %>% slice(grep("Steroid", pathway))
pathways$triacylglycerol <- ora_pathways %>% slice(grep("triacylglycerol", pathway))

pathways_keggIds <- pathways %>% 
  lapply(function(x) x %>% select(members_input_overlap) %>% unlist() %>% strsplit("; ") %>% 
           unlist() %>% unique())

intersect(pathways_keggIds$fatty_acids, pathways_keggIds$linolenic) # 8 ids overlap out of 18, and 9 ids
intersect(pathways_keggIds$fatty_acids, pathways_keggIds$phospholipases) # 11 ids overlap out of 18, and 13 ids

intersect(pathways_keggIds$amino_acids, pathways_keggIds$valine) # 5 ids overlap out of 28, and 9 ids
intersect(pathways_keggIds$amino_acids, pathways_keggIds$valine) # 5 ids overlap out of 28, and 17 ids

# use each pathway separately 
pathways_keggIds <- split(ora_pathways, f = ora_pathways$pathway ) %>%
  lapply(function(x) x %>% select(members_input_overlap) %>% unlist() %>% strsplit("; ") %>% 
           unlist() %>% unique())

# check the formulas --------------------------
library("KEGGREST")
library(tidyverse)

# load data ----------------------------------
load("cohorts_dat.RData")
#load("selected_DAMs.RData") # 45 formulas
load("selected_DAMs_padj2015.RData") # 146 formulas
# convert Formula to KEGG id ------------------
kegg_annot <- c() # KEGG database: using formula to annotate metabolite 
for (formula in selected_DAs) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
} # 45 formulas --> 126 Kegg IDs, 146 formulas --> 484 kegg IDs

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  separate(col = "cpdId", into = c("cpd", "keggIds"), sep = ":") 

# link Kegg id and formula -------------
pathways_formulas <- pathways_keggIds %>% 
  lapply(function(x) x %>% as.data.frame %>% 
           rename("keggIds" = ".") %>% left_join(kegg_outcomes) %>% 
           select(Formula) %>% unlist() %>% unique)

pathways_formulas_tab <- pathways_formulas %>% 
  lapply(function(x) x %>% as.data.frame %>% rename("formula" = ".")) %>%
  bind_rows(.id = "pathway")

# boxplot of the pathway expression level (avg. expression of all metabolite in the pathways)
## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

inputDat_pathway <- pathways_formulas_tab %>% 
  left_join(inputDat %>% column_to_rownames("name") %>% 
              t() %>% as.data.frame %>% rownames_to_column("formula")) %>%
  mutate(across(c(3:303), as.numeric)) %>% as.tibble()

inputDat_pathway_avg <- inputDat_pathway %>% select(-formula) %>% 
  group_by(pathway) %>% summarise(across(where(is.numeric), mean)) %>%
  column_to_rownames("pathway") %>% t() %>% as.data.frame %>% 
  rownames_to_column("name") %>% right_join(metadata_healthy)


# boxplot ================================
library(ggpubr)
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abT1 <- list(c("low", "high"))

mebo_pathway <- "fatty_acids"
mebo_pathway <- "linolenic" 
mebo_pathway <- "phospholipases"
mebo_pathway <- "amino_acids" 
mebo_pathway <-  "valine"
mebo_pathway <- "arginine"
mebo_pathway <- "sphingo" 
mebo_pathway <- "gAlpha" 
mebo_pathway <- "retinol"
mebo_pathway <- "steroid"
mebo_pathway <- "triacylglycerol"

mebo_pathway <- names(pathways_keggIds)[15]

metadat_boxplot <- inputDat_pathway_avg %>% 
  select(season, responder, c(mebo_pathway), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo_pathway,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

# based on abFC: NR vs R
ggboxplot(metadat_boxplot, x = "H1N1_abFC", y = mebo_pathway,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abFC, method = "t.test")

# based on abT1: low Ab vs high Ab
ggboxplot(metadat_boxplot, x = "H1N1_abBaseline", y = mebo_pathway,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abT1, method = "t.test")

# based on responders group: NR, other, TR
ggboxplot(metadat_boxplot, x = "responder", y = mebo_pathway,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_responder, method = "t.test")

