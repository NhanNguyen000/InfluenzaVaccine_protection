library(limma)
library(gplots)

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$name]
  
  if (identical(metaDat$name, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ + sex + age + reclassify, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

get.limmaRes_multiComparison <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with multiple comparison 
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limma model compare 2 groups (or 4 groups for LL_vsOther and HL_vsOther) together
  # res_4groups - list of outcome from limma model compare 4 groups separately
  
  res <- list()
  # res_4groups <- list()
  for (strain_group in strain_groups) {
    # LL vs LH, HL, HH together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na()
    
    res[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # # LL vs LH, HL, HH per group
    # reclass <- c("LH", "HL", "HH")
    # for (i in reclass) {
    #   metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #     select(name, sex, age, reclassify) %>% drop_na() %>%
    #     filter(reclassify %in% c("LL", i)) %>%
    #     mutate(reclassify = factor(reclassify, levels = c("LL", i)))
    #   
    #   res_4groups[[strain_group]]$LL_vsOthers[[i]] <- get.limmaRes(metadat_temp, inputDat)
    # }
    
    # LL vs (LH, HL, HH) as Protector
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = ifelse(reclassify == "LL", "LL", "Protector")) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", "Protector")))
    
    res[[strain_group]]$LL_vsProtectors <- get.limmaRes(metadat_temp, inputDat)
    
    # # LL vs (LH an HH) together
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>%
    #   filter(reclassify %in% c("LL", "LH", "HH")) %>%
    #   mutate(reclassify = ifelse(reclassify == "LL", "LL", "LH_HH")) %>%
    #   mutate(reclassify = factor(reclassify, levels = c("LL", "LH_HH")))
    # 
    # res[[strain_group]]$LL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)

    # # HL vs LH, HL, HH together
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>% 
    #   mutate(reclassify = factor(reclassify, levels = c("HL", "LL","LH", "HH")))
    # 
    # res[[strain_group]]$HL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # # HL vs (LH an HH) together
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>%
    #   filter(reclassify %in% c("HL", "LH", "HH")) %>%
    #   mutate(reclassify = ifelse(reclassify == "HL", "HL", "LH_HH")) %>%
    #   mutate(reclassify = factor(reclassify, levels = c("HL", "LH_HH")))
    # 
    # res[[strain_group]]$HL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
    
    # # HL vs HH 
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>%
    #   filter(reclassify %in% c("HL", "HH")) %>%
    #   mutate(reclassify = ifelse(reclassify == "HL", "HL", "HH"))  %>%
    #   mutate(reclassify = factor(reclassify, levels = c("HL", "HH")))
    
    # res[[strain_group]]$HL_vsHH <- get.limmaRes(metadat_temp, inputDat)
    
    # # (LL and HL as low reponse group) vs. (LH and HH as high response group)
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>%
    #   mutate(reclassify = substr(reclassify, 2, 2)) %>%
    #   mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    # res[[strain_group]]$LvsH_response <- get.limmaRes(metadat_temp, inputDat)

    # # (LL and LH as low baseline group) vs. (HL and HH as high baseline group)
    # metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
    #   select(name, sex, age, reclassify) %>% drop_na() %>%
    #   mutate(reclassify = substr(reclassify, 1, 1)) %>%
    #   mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    # res[[strain_group]]$LvsH_baseline <- get.limmaRes(metadat_temp, inputDat)
    
  }
  # return(list("res" = res, "res_4groups" = res_4groups))
  return(res)
}

# get.DE_tstat_vsOthers <- function(res) {
#   # Aim: extract the DE proteins/metabolites and the t-statistic value from 4 groups comparison ( LL_vsOther and HL_vsOther)
#   # following the get.limmaRes_multiComparesion's result
#   
#   resDE <- list()
#   res_tstatistic <- list()
#   for (compareType in c("LL_vsOthers", "HL_vsOthers")) {
#     compareGroups <- colnames(res[[1]][[compareType]]$p.value)[4:6]
#     
#     for (compareGroup in compareGroups) {
#       resDE[[compareType]][[compareGroup]] <- res %>%
#         lapply(function(x) x[[compareType]]$p.value %>% 
#                  as.data.frame() %>% select(compareGroup) %>% 
#                  filter(. <0.05))
#     }
#     
#     res_tstatistic[[compareType]] <- res %>%
#       lapply(function(x) x[[compareType]]$t %>% 
#                as.data.frame() %>% select(matches("reclassify")))
#   }
#   return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
# }

get.DE_tstat_LLvsOthers <- function(res) {
  # Aim: extract the DE proteins/metabolites and the t-statistic value from 4 groups comparison (LL_vsOther)
  # following the get.limmaRes_multiCompariison's result
  
  resDE <- list()
  res_tstatistic <- list()
  
  compareGroups <- colnames(res[[1]]$LL_vsOthers$p.value)[4:6]
  
  for (compareGroup in compareGroups) {
    resDE[[compareGroup]] <- res %>%
      lapply(function(x) x$LL_vsOthers$p.value %>%
               as.data.frame() %>% select(compareGroup) %>%
               filter(. <0.05))
  }
  
  res_tstatistic<- res %>%
    lapply(function(x) x$LL_vsOthers$t %>%
             as.data.frame() %>% select(matches("reclassify")))
  
  return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
}


# get.DE_tstat_vs1group <- function(res) {
#   # Aim: extract the DE proteins/metabolites and the t-statistic value from 2 groups comparison (not LL_vsOther, HL_vsOther)
#   # following the get.limmaRes_multiComparesion's result (res)
#   
#   resDE <- res %>%
#     lapply(function(x) x %>% 
#              purrr::list_modify("LL_vsOthers" = NULL, "HL_vsOthers" = NULL)) %>%
#     lapply(function(x) x %>% 
#              lapply(function(y) y$p.value %>% as.data.frame %>% 
#                       select(4) %>% filter(. < 0.05)))
#   
#   res_tstatistic <-  res %>%
#     lapply(function(x) x %>% 
#              purrr::list_modify("LL_vsOthers" = NULL, "HL_vsOthers" = NULL)) %>%
#     lapply(function(x) x %>% 
#              lapply(function(y) y$t %>% as.data.frame %>% select(4)) %>%
#              imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
#              lapply(function(y) y %>% rownames_to_column("valName")) %>%
#              purrr::reduce(full_join) %>% column_to_rownames("valName")) %>%
#     imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
#     lapply(function(x) x %>%
#              as.data.frame %>% rownames_to_column("valName")) %>% 
#     purrr::reduce(full_join)
#   
#   return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
# }

get.DE_tstat_LLvsProtector <- function(res) {
  # Aim: extract the DE proteins/metabolites and the t-statistic value from 2 groups comparison (LL vs. Protector)
  # following the get.limmaRes_multiComparison's result (res)
  
  resDE <- res %>%
    lapply(function(x) x$LL_vsProtectors$p.value %>% 
             as.data.frame %>% select(4) %>% filter(. < 0.05))
  
  res_tstatistic <-  res %>%
    lapply(function(x) x$LL_vsProtectors$t %>% 
             as.data.frame %>% select(4)) %>%
    imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
    lapply(function(y) y %>% rownames_to_column("valName")) %>%
    purrr::reduce(full_join) %>% column_to_rownames("valName")
  return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
}

# Protein data (use iMED is a discovery cohort, and ZirFlu as a replicate)-------------------------
# load data from 2 cohorts
load("../ZirFlu_NhanNguyen/ZirFlu.RData")
load("../iMED_NhanNguyen/iMED.RData")

## metadata for all healthy subjects -------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

## overlap proteins between iMED and ZirFlu -------------------------
overlapped_proteins <- intersect(iMED$protein_annot$OlinkID, ZirFlu$protein_annot$OlinkID)

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  filter(OlinkID %in% overlapped_proteins) %>%
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay")

## prepare data per cohort -------------------------
metadat_iMED <- metadat_healthy %>% filter(cohort == "iMED")
inputDat_iMED <- inputDat[, metadat_iMED$name]

metadat_ZirFlu_healthy <- metadat_healthy %>% filter(cohort == "ZirFlu")
inputDat_ZirFlu_healthy <- inputDat[, metadat_ZirFlu_healthy$name]

## run limma model for LL vs. protector, and LL vs. LH, HL, HH --------------------------
# iMED
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_iMED <- get.limmaRes_multiComparison(metadat =  metadat_iMED, 
                                          inputDat = inputDat_iMED, 
                                          strain_groups = strain_groups)

# ZirFlu
res_ZirFlu_healthy <- get.limmaRes_multiComparison(metadat =  metadat_ZirFlu_healthy, 
                                                    inputDat = inputDat_ZirFlu_healthy, 
                                                    strain_groups = "H1N1_reclassify")

# Metabolite data (use iMED is a discovery cohort, and ZirFlu as a replicate)-------------------------
## metadata for all healthy subjects -------------------------
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

## iMED metabolites -------------------------
iMED_metabolites <- iMED$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

iMED_metaboliteDat <- iMED$metabolite
if (identical(iMED_metabolites$ionIdx, as.numeric(colnames(iMED$metabolite))) == TRUE) {
  colnames(iMED_metaboliteDat) <- iMED_metabolites$Formula
}


## ZirFlu metabolites -------------------------
ZirFlu_metabolites <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

ZirFlu_metaboliteDat <- ZirFlu$metabolite_dat
if (identical(ZirFlu_metabolites$ionIdx, as.numeric(colnames(ZirFlu$metabolite_dat))) == TRUE) {
  colnames(ZirFlu_metaboliteDat) <- ZirFlu_metabolites$Formula
}

## overlap metabolites between iMED and ZirFlu -------------------------
overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

metadat_iMED <- metadat_healthy %>% filter(cohort == "iMED")
inputDat_iMED <- iMED_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_iMED$name)

metadat_ZirFlu_healthy <- metadat_healthy %>% filter(cohort == "ZirFlu")
inputDat_ZirFlu_healthy <- ZirFlu_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_ZirFlu_healthy$name)

## run limma model for LL vs. protector (LH, HL, HH), and LL vs. LH --------------------------
# iMED
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_iMED <- get.limmaRes_multiComparison(metadat =  metadat_iMED, 
                                          inputDat = inputDat_iMED, 
                                          strain_groups = strain_groups)

# ZirFlu
res_ZirFlu_healthy <- get.limmaRes_multiComparison(metadat =  metadat_ZirFlu_healthy, 
                                                    inputDat = inputDat_ZirFlu_healthy, 
                                                    strain_groups = "H1N1_reclassify")
