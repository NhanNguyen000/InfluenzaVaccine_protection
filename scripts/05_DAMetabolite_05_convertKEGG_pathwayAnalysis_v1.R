rm(list = ls())

library("KEGGREST")
library(tidyverse)

# load data ----------------------------------
load("cohorts_dat.RData")
load("selected_DAMs.RData")

# convert Formula to KEGG id ------------------
## selected DAMs -------------------
kegg_annot <- c() # KEGG database: using formula to annotate metabolite 
for (formula in selected_DAs) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
} # 45 formulas --> 126 Kegg IDs

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  separate(col = "cpdId", into = c("cpd", "keggIds"), sep = ":") 
length(unique(kegg_outcomes$Formula)) # 38 formulas were able to annotate

write.table(kegg_outcomes$keggIds, 
            file = "keggID_DAMs_reclassify_4to2groups.txt",
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
