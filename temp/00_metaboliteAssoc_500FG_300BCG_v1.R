# Martijn data: /vol/projects/CIIM/Influenza/iMED/supplemental_data  
# Associations of 300BCG and 500FG metabolites to cytokines and cellcounts in the subfolders here. Look in output/assoc_met_*.csv

C4H8O3_annot <- iMED$metabolite_annot %>% filter(Formula == "C4H8O3")

# check in 500FG
assocTypes <- c("assoc_met_cellcount.csv", "assoc_met_cytokine.csv", "assoc_met_prot.csv")
for (assocType in assocTypes) {
  dat <- read.csv(paste0("/vol/projects/CIIM/Influenza/iMED/supplemental_data/500FG/output/", 
                         assocType))
  C4H8O3_assoc_500FG[[assocType]] <- dat %>% 
    filter(Top.annotation.ids %in% C4H8O3_annot$CompoundID) %>% 
    filter(p.value < 0.05)
} # no C4H8O3 in 500FG

# check in 300BCG
C4H8O3_assoc_300BCG <- list()
assocTypes <- c("assoc_met_cellcounts.csv", "assoc_met_cytokines.csv", "assoc_met_prot.csv")
for (assocType in assocTypes) {
  dat <- read.csv(paste0("/vol/projects/CIIM/Influenza/iMED/supplemental_data/300BCG/output/", 
                         assocType))
  
  C4H8O3_assoc_300BCG[[assocType]] <- dat %>% 
    filter(Top.annotation.ids %in% C4H8O3_annot$CompoundID) %>% 
    filter(p.value < 0.05)
}


