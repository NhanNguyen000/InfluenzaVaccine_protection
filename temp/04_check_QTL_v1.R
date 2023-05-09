# eQTL ---------------------------------
# from Hua
dat <- read.table("data/CD83_single_cell_eqtl.txt", header = TRUE)
dat_pvalue <- dat %>% filter(p_value < 0.05)

CD83_B <- dat_pvalue %>% filter(cell_type == "B")

# form Javi
library(phenoscanner)

df_list <- split(dat, rep(1:ceiling(nrow(dat)/100), each = 100, length.out = nrow(dat)))
outcome <- list()

for (i in names(df_list)) {
  loci <- df_list[[i]] %>% dplyr::rename('rsid' = snp_id, 'pval'= p_value)
  
  #Add phenoscanner info about gwas, eqtl and annotation
  ps <- phenoscanner(loci$rsid, catalogue = "GWAS")
  pos <- ps$snps %>%
    select(c("snp","consequence","hgnc"))%>%
    rename("rsid" = snp)
  eqtl <- phenoscanner(loci$rsid, catalogue = "eQTL")$results %>%
    select(c("snp","exp_gene"))%>%
    rename("rsid" = snp)
  gwas <- ps$results %>%
    select(c("snp","trait"))%>%
    rename("rsid" = snp)
  
  outcome[[i]] <- loci %>%
    left_join(pos, "rsid")%>%
    left_join(eqtl, "rsid")%>%
    group_by(across(c(-exp_gene)))%>%
    summarise_at("exp_gene", function(x) paste(unique(x), collapse = ";"))%>%
    left_join(gwas, "rsid")%>%
    group_by(across(c(-trait)))%>%
    summarise_at("trait", function(x) paste(unique(x), collapse = ";"))
  
}

loci_info <- outcome %>% purrr::reduce(full_join)
unique(loci_info$trait) # not things related to Influenza
unique(loci_info$hgnc) # not things related to Influenza

# pQTL from Matijn --------------------------------
## pQTL from BioBank -------------------
# /vol/projects/CIIM/resources/pQTL_biobank/meta_filtered_final, directory
protein_withQTL <- list.dirs("/vol/projects/CIIM/resources/pQTL_biobank/meta_filtered_final")
grep("CD83", protein_withQTL) # no pQTL of CD83
grep("BL11", protein_withQTL) # other name of CD83
grep("HB15", protein_withQTL) # other name of CD83 
protein_withQTL[grep("/CD", protein_withQTL)] # have pQTL of CD38, CD84, CD93, ...

## pQTL form iMED -----------------------
pQTL_iMED_CD83 <- read.table("/vol/projects/CIIM/Influenza/iMED/proteomic/pQTL/time_individual/output/T1/outs/out_phenotype_prot_89.csv",
                        header = TRUE)
# gene == "OID20565",  OlinkID of CD83

range(pQTL_iMED_CD83$p.value)
pQTL_iMED_CD83_pvalue <- pQTL_iMED_CD83 %>% filter(p.value < 0.05)
