# load libraries & functions --------------------------
library(fgsea)
#library('org.Hs.eg.db')
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

get.fgsea_singleInput <- function(dat, col_input, col_geneID) {
  # Description: prepare the input for fgsea
  
  # Arguments: 
  # dat - data with gene ID in one column and gene value/expression in other columns
  # col_input - name of the input column
  # col_geneID - name of the column present gene ID
  
  # Returns: fgsaInput - input vector for the fgsea analysis
  
  fgseaInput <- dat[, col_input]
  names(fgseaInput) <- dat[, col_geneID]
  
  return(fgseaInput)
}

get.fgseaRes <- function(dat, convertID, hallmarkSet) {
  # Description: run fgsea (GSEA) for all subjects/patients
  
  # Arguments: 
  # dat - data with patients in row, proteins in column
  # convertID - table to convert protein OlinkID to gene Entrez ID
  # hallmarkSet - the hallmark / reference gene set to run GSEA
  
  # Returns: fgsaRes - a list of fgsea result for each patient
  
  fgseaRes <- list()
  for (patient in rownames(dat)) {
    
    ranks_temp <- get.fgsea_singleInput(
      dat = dat %>% t() %>% as.data.frame() %>% 
        rownames_to_column("OlinkID") %>% left_join(convertID), 
      col_input = patient, col_geneID = "ENTREZID")
    
    fgseaRes[[patient]] <- fgsea(hallmarkSet, ranks_temp) %>% 
      as_tibble() %>% dplyr::select(pathway, size, leadingEdge, NES) %>%
      mutate(name = patient)
  }
  
  return(fgseaRes)
}


# calculate inflammation score - using average (NES values of inflammation hallmark sets from GSEA) at baseline  -----------------------
load("proteinDat_impute.RData") # load imputed protein data
load("gsea_hallmarktSet.RData") # load hallmark gene sets (prepare from from MSigDB) & convert protein - Entrez gene ID

# run code for all cohort
fgseaResTab <- list()
for (cohort in names(proteinDat_impute)) {
  inputDat <- proteinDat_impute[[cohort]]
  
  fgseaRes <- list()
  fgseaRes[["hallmarkSets"]] <- get.fgseaRes(dat = inputDat,
                                             convertID = protein_Entrez, 
                                             hallmarkSet = hallmarkSets) 
  
  fgseaRes[["hallmarkSubsets_v1"]] <- get.fgseaRes(dat = inputDat,
                                                   convertID = protein_Entrez,
                                                   hallmarkSet = hallmarkSets[hallmarkSubsets_v1]) 
  
  fgseaRes[["hallmarkSubsets_v2"]] <- get.fgseaRes(dat = inputDat,
                                                   convertID = protein_Entrez,
                                                   hallmarkSet = hallmarkSets[hallmarkSubsets_v2]) 
  
  fgseaResTab[[cohort]] <- fgseaRes %>% 
    lapply(function(x) x %>% purrr::reduce(full_join) %>% # long format
             dplyr::select(-size, -leadingEdge) %>%
             pivot_wider(names_from = name, values_from = NES) %>% # convert to wide format
             column_to_rownames("pathway"))
  
}

inflamScore <- list()
for (cohort in names(fgseaResTab)) {
  inflamScore[[cohort]] <- fgseaResTab[[cohort]] %>% 
    lapply(function(x) x %>% t() %>% as.data.frame() %>%
             dplyr::select(all_of(hallmarkSubsets_v2)) %>%  # inflamScore use hallmarkSubsets_v2 pathways
             rowMeans() ) %>%
    as.data.frame() %>% rename_with(~paste0("inflamScore.", .x))
}

# Checked: Not much different in the NES result when run all hallmarset, 4/6 inflammatory sets 

# check if things are correct  ------------------------------
identical(metadat_sampleT1$ZirFlu_2019$probenID, colnames(fgseaResTab$ZirFlu_2019$hallmarkSets)) # TRUE
identical(rownames(inflamScore$ZirFlu_2019), names(fgseaResTab$ZirFlu_2019$hallmarkSets)) # TRUE

# inflammation score across age ----------------------------------------------------------------------
inflamScore.NESfgsea <- inflamScore %>% 
  purrr::reduce(rbind) %>% rownames_to_column("name") %>%
  left_join(cohorts_dat$donorSample_all)

save(fgseaResTab, inflamScore.NESfgsea, file = "inflamScore_NESfgsae.RData")
