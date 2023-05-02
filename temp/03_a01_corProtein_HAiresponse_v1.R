# check correlation between protein vs. HAI response -----
library(rstatix)
ZirFlu_metadatT1 <- ZirFlu$HAItiter %>% filter(season == "2019") %>% 
  left_join(ZirFlu$donorSamples %>% filter(season == "2019" & time == "T1"))

cor(ZirFlu_proteinDat$OID20429, ZirFlu_metadatT1$H1N1_T1, )

k <- ZirFlu_metadatT1 %>% full_join(ZirFlu_proteinDat %>% rownames_to_column("probenID"))
cor_res <- cor_test(k %>% dplyr::select(OID20429, H1N1_T1))

ZirFlu_metadatT1_v2 <- ZirFlu_metadatT1 %>% 
  select(probenID, matches("T1|T2|T3|abFC")) %>% column_to_rownames("probenID")

outcome <- c()
for (val in names(ZirFlu_metadatT1_v2)) {
  outcome_temp <- c()
  
  for (protein in names(ZirFlu_proteinDat)) {
    dat_temp <- ZirFlu_metadatT1_v2 %>% 
      dplyr::select(any_of(val)) %>% rownames_to_column("matched_keys") %>%
      full_join(ZirFlu_proteinDat %>% 
                  dplyr::select(any_of(protein)) %>% rownames_to_column("matched_keys")) %>%
      column_to_rownames("matched_keys")
    
    cor_res <- cor_test(dat_temp)
    outcome_temp <- rbind(outcome_temp, cor_res)
  }
  outcome[[val]] <- outcome_temp
}
View(outcome$H1N1_T1)


fsea_datInput <- outcome$H1N1_T1 %>% select(var2, cor) %>% 
  left_join(protein_Entrez, by = c("var2" = "OlinkID"))

ranks <- fsea_datInput$cor
names(ranks) <- fsea_datInput$ENTREZID
#barplot(sort(ranks, decreasing = T))
