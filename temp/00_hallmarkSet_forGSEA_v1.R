# load libraries & functions --------------------------
library(fgsea)
library('org.Hs.eg.db')

# load hallmark gene sets from MSigDB --------------------------
hallmarkSets <- gmtPathways("reference/h.all.v2022.1.Hs.entrez.gmt.txt")
# BTM_sets <- gmtPathways("reference/BTM_for_GSEA_20131008.gmt")
hallmarkSets %>% head() %>% lapply(head)

# selected hallmark sets
all_hallmarkSets <- names(hallmarkSets)

hallmarkSubsets_v1 <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING", 
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_COMPLEMENT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

hallmarkSubsets_v2 <- c(
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_COMPLEMENT", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

# convert protein to Entrez gene ID
Entrez_IDs <- select(org.Hs.eg.db, 
                     keys = cohorts_dat$proteinAnnot_all$Assay,
                     keytype = 'SYMBOL', columns = c('ENTREZID', 'SYMBOL'))

protein_Entrez <- cohorts_dat$proteinAnnot_all %>% 
  full_join(Entrez_IDs, by = c("Assay" = "SYMBOL"))

protein_Entrez$ENTREZID[protein_Entrez$Assay == "CXCL14"] <- 9547
protein_Entrez$ENTREZID[protein_Entrez$Assay == "ALDH3A1"] <- 218
protein_Entrez$ENTREZID[protein_Entrez$Assay == "TRIM5"] <- 85363
protein_Entrez$ENTREZID[protein_Entrez$Assay == "MICB_MICA"] <- 100507436 # or 4277 both of them are not in the inflammatory hallmark sets
protein_Entrez$ENTREZID[protein_Entrez$Assay == "CKMT1A_CKMT1B"] <- 548596 # or 1159 both of them are not in the inflammatory hallmark sets

# save data -------------------------------------------
save(hallmarkSets, hallmarkSubsets_v1, hallmarkSubsets_v2, protein_Entrez, file = "gsea_hallmarktSet.RData")

