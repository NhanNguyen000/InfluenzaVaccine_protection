# run thi code in the terminal ---------------------------
# because the vcf file have some prefix in the SNPs, so could not use the bcftools, we use the defaut linux command
cd /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/
grep -w -e CHROM -f CXCL1_rsid.txt input/iMED_vcf.vcf >| CXCL1_genotype.txt 

# run the rest of the code in R ---------------------------
SNPs <- read.table("processedDat/omicsPred/CXCL1_rsid.txt") %>% unlist()
genotype_raw <- read.table("processedDat/omicsPred/CXCL1_genotype.txt", header = TRUE)

genotype <- genotype_raw %>% 
  mutate(across(c(10:168), function(x) paste0(substring(x, 1, 1), substring(x, 3, 3)))) %>%
  mutate(across(c(10:168), 
                function(x) ifelse(x == "00", paste0(REF, REF),
                                   ifelse(x == "11", paste0(ALT, ALT),
                                          ifelse(x == "10" |x == "01", paste0(REF, ALT),
                                                 ifelse(x == "..", "..", NA)))))) %>%
  separate(col = ID, into = c("CHR:POS", "rsid"), sep = "%") %>% 
  select(rsid, matches("I.")) %>% select(-2, -3) %>%
  column_to_rownames("rsid") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("probandID") %>% mutate(probandID = gsub("[.]", "-", probandID))

# use normalized(abTiter T4) - Martijn did 
normalized_abTiter <- read.table("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/phenotypes.csv") %>% 
  t() %>% as.data.frame %>% rownames_to_column("probandID") %>%
  mutate(probandID = gsub("X", "", probandID)) %>%
  mutate(probandID = ifelse(nchar(probandID) == 1, 
                            paste0("I-000", probandID), 
                            ifelse(nchar(probandID) == 2, paste0("I-00", probandID), 
                                   ifelse(nchar(probandID) == 3, 
                                          paste0("I-0", probandID), paste0("I-", probandID)))))

genotype_abTiter <- genotype %>% left_join(normalized_abTiter)

plotList <- list()
for (rsid in SNPs) {
  plotList[[paste0(rsid, "_H1N1")]] <- genotype_abTiter %>% 
    ggplot(aes_string(x = rsid, y = "H1N1_T4")) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
    theme_bw()
  plotList[[paste0(rsid, "_H3N2")]] <- genotype_abTiter %>% 
    ggplot(aes_string(x = rsid, y = "H3N2_T4")) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
    theme_bw()
  plotList[[paste0(rsid, "_B")]] <- genotype_abTiter %>% 
    ggplot(aes_string(x = rsid, y = "B_T4")) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
    theme_bw()
}

gridExtra::grid.arrange(grobs = plotList, ncol = 6, as.table = FALSE)

