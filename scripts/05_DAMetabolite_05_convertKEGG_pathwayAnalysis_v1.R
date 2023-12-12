rm(list = ls())

library("KEGGREST")
library(tidyverse)

# load data ----------------------------------
load("cohorts_dat.RData")
#load("selected_DAMs.RData") # 45 formulas
#load("selected_DAMs_padj2015.RData") # 146 formulas
load("selected_consist_sigPval_DAMs.RData") # 30 and 180 formulas
load("sigMebo_dynamic.RData")

# convert Formula to KEGG id ------------------
## selected DAMs -------------------
kegg_annot <- c() # KEGG database: using formula to annotate metabolite 
selected_DAs <- selected_vars
selected_DAs <- sigPval_vars
selected_DAs <- unique(sigOutcome$var)

for (formula in selected_DAs) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
} # 45 formulas --> 126 Kegg IDs, 146 formulas --> 484 kegg IDs, 70 fomulas --> 223 kegg IDs

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  separate(col = "cpdId", into = c("cpd", "keggIds"), sep = ":") 
length(unique(kegg_outcomes$Formula)) # 38 / 116 / 58 / 25 / 145 formulas were able to annotate

write.table(kegg_outcomes$keggIds, 
            #file = "keggID_DAMs_reclassify_4to2groups.txt",
            #file = "keggID_DAMs_reclassify_4groups.txt",
            #file = "keggID_consist_sigPval_DAMs.txt", # 25 formulas to 107 kegg ID  were able to annotate
            #file = "keggID_all_sigPval_DAMs.txt", # 145 formulas to 670 kegg ID were able to annotate
            file = "keggID_all_sigMeboDaynamic.txt", # 265 formulas to 1198 kegg ID were able to annotate
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## all the endogenous metabolites used in the lm() model --------------------------
kegg_annot <- c() # KEGG database: using formula to annotate metabolite

all_Formulas <- colnames(mebo_Dat$iMED_2014)

for (formula in unique(all_Formulas)) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
} # 508 formulas --> 1731 Ids

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  separate(col = "cpdId", into = c("cpd", "keggIds"), sep = ":") 
length(unique(kegg_outcomes$Formula)) # 376 formulas were able to annotate

write.table(kegg_outcomes$keggIds, 
            file = "keggID_allMetabolites.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
