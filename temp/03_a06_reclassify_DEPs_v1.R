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

get.limmaRes_multiComparesion <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with multiple comparison 
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limmar model compare 2 groups (or 4 groups for LL_vsOther and HL_vsOther) together
  # res_4groups - list of outcome from limmar model compare 4 groups separately
  
  res <- list()
  res_4groups <- list()
  for (strain_group in strain_groups) {
    # LL vs LH, HL, HH together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na()
    
    res[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # LL vs LH, HL, HH per group
    reclass <- c("LH", "HL", "HH")
    for (i in reclass) {
      metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
        select(name, sex, age, reclassify) %>% drop_na() %>%
        filter(reclassify %in% c("LL", i)) %>%
        mutate(reclassify = factor(reclassify, levels = c("LL", i)))
      
      res_4groups[[strain_group]]$LL_vsOthers[[i]] <- get.limmaRes(metadat_temp, inputDat)
    }
    
    # LL vs (LH, HL, HH) as Protector
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = ifelse(reclassify == "LL", "LL", "Protector")) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", "Protector")))
    
    res[[strain_group]]$LL_vsProtectors <- get.limmaRes(metadat_temp, inputDat)
    
    # LL vs (LH an HH) together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("LL", "LH", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "LL", "LL", "LH_HH")) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", "LH_HH")))
    
    res[[strain_group]]$LL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs LH, HL, HH together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>% 
      mutate(reclassify = factor(reclassify, levels = c("HL", "LL","LH", "HH")))
    
    res[[strain_group]]$HL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs (LH an HH) together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("HL", "LH", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "HL", "HL", "LH_HH")) %>%
      mutate(reclassify = factor(reclassify, levels = c("HL", "LH_HH")))
    
    res[[strain_group]]$HL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs HH 
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("HL", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "HL", "HL", "HH"))  %>%
      mutate(reclassify = factor(reclassify, levels = c("HL", "HH")))
    
    res[[strain_group]]$HL_vsHH <- get.limmaRes(metadat_temp, inputDat)
    
    # (LL and HL as low reponse group) vs. (LH and HH as high response group)
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = substr(reclassify, 2, 2)) %>%
      mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    res[[strain_group]]$LvsH_response <- get.limmaRes(metadat_temp, inputDat)
    
    # (LL and LH as low baseline group) vs. (HL and HH as high baseline group)
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = substr(reclassify, 1, 1)) %>%
      mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    res[[strain_group]]$LvsH_baseline <- get.limmaRes(metadat_temp, inputDat)
    
  }
  return(list("res" = res, "res_4groups" = res_4groups))
}

get.DE_tstat_vsOthers <- function(res) {
  # Aim: extract the DE proteins/metabolites and the t-statistic value from 4 groups comparison ( LL_vsOther and HL_vsOther)
  # following the get.limmaRes_multiComparesion's result
  
  resDE <- list()
  res_tstatistic <- list()
  for (compareType in c("LL_vsOthers", "HL_vsOthers")) {
    compareGroups <- colnames(res[[1]][[compareType]]$p.value)[4:6]
    
    for (compareGroup in compareGroups) {
      resDE[[compareType]][[compareGroup]] <- res %>%
        lapply(function(x) x[[compareType]]$p.value %>% 
                 as.data.frame() %>% select(compareGroup) %>% 
                 filter(. <0.05))
    }
    
    res_tstatistic[[compareType]] <- res %>%
      lapply(function(x) x[[compareType]]$t %>% 
               as.data.frame() %>% select(matches("reclassify")))
  }
  return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
}

get.DE_tstat_vs1group <- function(res) {
  # Aim: extract the DE proteins/metabolites and the t-statistic value from 2 groups comparison (not LL_vsOther, HL_vsOther)
  # following the get.limmaRes_multiComparesion's result (res)
  
  resDE <- res %>%
    lapply(function(x) x %>% 
             purrr::list_modify("LL_vsOthers" = NULL, "HL_vsOthers" = NULL)) %>%
    lapply(function(x) x %>% 
             lapply(function(y) y$p.value %>% as.data.frame %>% 
                      select(4) %>% filter(. < 0.05)))
  
  res_tstatistic <-  res %>%
    lapply(function(x) x %>% 
             purrr::list_modify("LL_vsOthers" = NULL, "HL_vsOthers" = NULL)) %>%
    lapply(function(x) x %>% 
             lapply(function(y) y$t %>% as.data.frame %>% select(4)) %>%
             imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
             lapply(function(y) y %>% rownames_to_column("valName")) %>%
             purrr::reduce(full_join) %>% column_to_rownames("valName")) %>%
    imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
    lapply(function(x) x %>%
             as.data.frame %>% rownames_to_column("valName")) %>% 
    purrr::reduce(full_join)
  
  return(list("resDE" = resDE, "res_tstatistic" = res_tstatistic))
}

get.plotDat_clusterRow <- function(plotDat, colName, valColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function witll write the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "valName" colum, and rowname from fixed "valName" column
  
  # outcome: plotDat_order - plot data long format with order in valName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% select(c("valName", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("valName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(valName = factor(valName, levels = valName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}

# input data (use iMED is a discovery cohort, and ZirFlu as a replicate)-------------------------
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

# run limma model for LL vs. protector (LH, HL, HH), and LL vs. LH --------------------------
# iMED
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_iMED <- get.limmaRes_multiComparesion(metadat =  metadat_iMED, 
                                          inputDat = inputDat_iMED, 
                                          strain_groups = strain_groups)

# ZirFlu
res_ZirFlu_healthy <- get.limmaRes_multiComparesion(metadat =  metadat_ZirFlu_healthy, 
                                                    inputDat = inputDat_ZirFlu_healthy, 
                                                    strain_groups = "H1N1_reclassify")

# sigProteins ---------------------------------------------------------
# LL or HL vs other groups comparison
k<- get.DE_tstat_vsOthers(res = res_iMED$res)
resDE_vsOthers <- k$resDE
resDE_vsOthers_tstatistic <- k$res_tstatistic

# LL or HL vs. other group in 1:1 group comparison
k2 <- get.DE_tstat_vs1group(res = res_iMED$res)
resDE_vs1group <- k2$resDE
resDE_vs1group_tstatistic <- k2$res_tstatistic

# LL vs Other per 1:1 group comparison
resSig_LLvsOther_perGroup <- res_iMED$res_4groups %>%
  lapply(function(x) x$LL_vsOther %>% 
           lapply(function(y) y$p.value %>% 
                    as.data.frame() %>% select(4) %>% filter(. < 0.05)))

resSig_LLvsOther_perGroup_tstatistic <- res_iMED$res_4groups %>%
  lapply(function(x) x$LL_vsOthers %>% 
           lapply(function(y) y$t %>% as.data.frame() %>% 
                    select(4) %>% rownames_to_column("Assay")) %>% 
           purrr::reduce(full_join) %>% relocate(Assay)
  )

# Check venn diagram -----------------------------------------------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = rownames(resDE_vs1group$H1N1_reclassify$LL_vsProtectors),
  "H3N2" = rownames(resDE_vs1group$H3N2_reclassify$LL_vsProtectors),
  "B" = rownames(resDE_vs1group$B_reclassify$LL_vsProtectors))
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
H1N1_B_proteins <- intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)

# heatmap  ----------------------------------------------------
# t-statistic data
tstat_vsOthers_wideDat <- resDE_vsOthers_tstatistic %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% as.data.frame %>% rownames_to_column("valName") %>%
                    gather(key = "compare", value = "tstatistic", -1) ))

tstat_vsOthers_longDat <- tstat_vsOthers_wideDat %>% 
  lapply(function(x) x %>% imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join))

tstat_vs1group_longDat <- resDE_vs1group_tstatistic %>% 
  gather(key = "compare", value = "tstatistic", -1) %>%
  separate(compare, sep = "[.]", into = c("strain", "compare", "compareName"))

# sig DEP between LL or HL vs. Others group
DE_vsOthers_wideDat <- resDE_vsOthers %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
                    lapply(function(z) z %>% as.data.frame %>% rownames_to_column("valName")) %>% 
                    purrr::reduce(full_join)) %>%
           purrr::reduce(full_join))

DE_vsOthers_longDat <- DE_vsOthers_wideDat %>%
  lapply(function(x) x %>% gather(key = "compare", value = "sig.p.value", -1) %>%
           separate(compare, sep = "[.]", into = c("strain", "compare")) %>%
           drop_na(sig.p.value))

DE_vs1group_longDat <- resDE_vs1group %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rename("sig.p.value" = 1) %>% rownames_to_column("valName")) %>%
           imap(~mutate(.x, compare = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, strain = .y)) %>%
  purrr::reduce(full_join)

# part 1: "LL_vsProtectors", "LvsH_response", "LvsH_baseline"  
compareType <- "LL_vsProtectors" #"LL_vsLH_HH", "HL_vsLH_HH", "HL_vsHH", "LvsH_response", "LvsH_baseline"  
# compareType <- "LvsH_response"
# compareType <- "LvsH_baseline"

DE_temp <- DE_vs1group_longDat %>% filter(compare == compareType)

plotDat <- tstat_vs1group_longDat %>% full_join(DE_vs1group_longDat) %>%
  filter(compare == compareType, valName %in% DE_temp$valName)

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "strain", 
                                        valColumn = "tstatistic")
plotDat_order %>%
  ggplot(aes(x = strain, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order %>%
  ggplot(aes(x = strain, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_flip()
# part 2: "LvsH_response"  

# next step:
# pick strain & plot heatmap
compareRef <- "LL_vsOthers" # "HL_vsOthers"
#strainType <- "H1N1_reclassify" # "H3N2_reclassify", "B_reclassify"

DE_temp <- DE_vsOthers_longDat[[compareRef]] %>% drop_na(sig.p.value)
DEs <- DE_temp$valName

DEs <- unique(c(intersect(res.venn$H1N1, res.venn$H3N2),
                intersect(res.venn$H1N1, res.venn$B),
                intersect(res.venn$H3N2, res.venn$B))) # DEPs in at least 2 strains
DEs <- res.venn$H1N1[-which(res.venn$H1N1 %in% c(res.venn$H3N2, res.venn$B))] # DEPs only in H1N1
DEs <- res.venn$B[-which(res.venn$B %in% c(res.venn$H3N2, res.venn$H1N1))] # DEPs only in B
DEs <- res.venn$H3N2[-which(res.venn$H3N2 %in% c(res.venn$B, res.venn$H1N1))] # DEPs only in H3N2

plotDat <- tstat_vsOthers_longDat[[compareRef]] %>% 
  full_join(DE_vsOthers_longDat[[compareRef]]) %>%
  filter(valName %in% DEs) %>%
  relocate(strain, .before = "compare") %>%
  unite("compare", strain:compare) %>% 
  mutate(compare = gsub("reclassify_reclassify", "LLvs", compare))


plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "compare", 
                                        valColumn = "tstatistic")

plotDat_order %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# LL vs Other (only 1 strain)
DEs_temp <- unique(unlist(res.venn)) # all DEs 
DEs_v2 <- DEs_temp[-which(DEs_temp %in% DEs)] # DEs in only 1 strain

plotDat <- tstat_vsOthers_longDat[[compareRef]] %>% 
  full_join(DE_vsOthers_longDat[[compareRef]]) %>%
  filter(valName %in% DEs_v2) %>%
  relocate(strain, .before = "compare") %>%
  unite("compare", strain:compare) %>% 
  mutate(compare = gsub("reclassify_reclassify", "LLvs", compare))


plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "compare", 
                                        valColumn = "tstatistic")

plotDat_order %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


plotDat_order_iMED <- plotDat_order %>%
  mutate(compare = paste0("iMED_", compare))
# validataion in ZirFlu ----------------
# LL or HL vs other groups comparison
g <- get.DE_tstat_vsOthers(res = res_ZirFlu_healthy$res)
ZirFlu_resDE_vsOthers <- g$resDE
ZirFlu_resDE_vsOthers_tstatistic <- g$res_tstatistic

# LL or HL vs. other group in 1:1 group comparison
g2 <- get.DE_tstat_vs1group(res = res_ZirFlu_healthy$res)
ZirFlu_resDE_vs1group <- g2$resDE
ZirFlu_resDE_vs1group_tstatistic <- g2$res_tstatistic

# t-statistic data
ZirFlu_tstat_vsOthers_wideDat <- ZirFlu_resDE_vsOthers_tstatistic %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% as.data.frame %>% rownames_to_column("valName") %>%
                    gather(key = "compare", value = "tstatistic", -1) ))

ZirFlu_tstat_vsOthers_longDat <- ZirFlu_tstat_vsOthers_wideDat %>% 
  lapply(function(x) x %>% imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join))

ZirFlu_tstat_vs1group_longDat <- ZirFlu_resDE_vs1group_tstatistic %>% 
  gather(key = "compare", value = "tstatistic", -1) %>%
  separate(compare, sep = "[.]", into = c("strain", "compare", "compareName"))

# part 1: "LL_vsProtectors", "LvsH_response", "LvsH_baseline"  
compareType <- "LL_vsProtectors" #"LL_vsLH_HH", "HL_vsLH_HH", "HL_vsHH", "LvsH_response", "LvsH_baseline"  
# compareType <- "LvsH_response"
# compareType <- "LvsH_baseline"

ZirFlu_DE_vsOthers_wideDat <- ZirFlu_resDE_vsOthers %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    imap(~.x %>% rename_with(function(x) paste(.y, x, sep = "."))) %>%
                    lapply(function(z) z %>% as.data.frame %>% rownames_to_column("valName")) %>% 
                    purrr::reduce(full_join)) %>%
           purrr::reduce(full_join))

ZirFlu_DE_vsOthers_longDat <- ZirFlu_DE_vsOthers_wideDat %>%
  lapply(function(x) x %>% gather(key = "compare", value = "sig.p.value", -1) %>%
           separate(compare, sep = "[.]", into = c("strain", "compare")) %>%
           drop_na(sig.p.value))

ZirFlu_DE_vs1group_longDat <- ZirFlu_resDE_vs1group %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rename("sig.p.value" = 1) %>% rownames_to_column("valName")) %>%
           imap(~mutate(.x, compare = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, strain = .y)) %>%
  purrr::reduce(full_join)

# ZirFLu plot dat
ZirFlu_plotDat <- ZirFlu_tstat_vs1group_longDat %>% 
  full_join(ZirFlu_DE_vs1group_longDat) %>%
  filter(compare == compareType, valName %in% levels(plotDat_order$valName)) %>%
  mutate(strain = "H1N1_ZirFlu")


plotDat_order_v2 <- plotDat_order %>%
  full_join(ZirFlu_plotDat) %>%
  mutate(valName = factor(valName, levels = levels(plotDat_order$valName), ordered = TRUE)) %>%
  mutate(strain = factor(strain, levels = c("B_reclassify", "H1N1_reclassify", "H3N2_reclassify", "H1N1_ZirFlu")))

plotDat_order_v2 %>%
  ggplot(aes(x = strain, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order_v3 <- get.plotDat_clusterRow(plotDat_order_v2, 
                                           colName = "strain", valColumn = "tstatistic")
plotDat_order_v3 %>%
  ggplot(aes(x = strain, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order_v3.2 <- plotDat_order_v3 %>% 
  mutate(strain = gsub("_reclassify", "_iMED", strain)) %>%
  mutate(strain = factor(strain, levels = c("H1N1_ZirFlu", "B_iMED", "H3N2_iMED", "H1N1_iMED")))

plotDat_order_v3.2 %>%
  ggplot(aes(x = valName, y = strain, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# part 2: LL vs. Other (LH, HL, HH)
# pick strain & plot heatmap
compareRef <- "LL_vsOthers" # "HL_vsOthers"
#strainType <- "H1N1_reclassify" # "H3N2_reclassify", "B_reclassify"

DEs <- unique(c(intersect(res.venn$H1N1, res.venn$H3N2),
                intersect(res.venn$H1N1, res.venn$B),
                intersect(res.venn$H3N2, res.venn$B))) # DEs in at least 2 strains

plotDat <- ZirFlu_tstat_vsOthers_longDat[[compareRef]] %>% 
  full_join(DE_vsOthers_longDat[[compareRef]]) %>%
  filter(valName %in% DEs) %>%
  relocate(strain, .before = "compare") %>%
  unite("compare", strain:compare) %>% 
  mutate(compare = gsub("reclassify_reclassify", "LLvs", compare))

ZirFlu_plotDat <- ZirFlu_tstat_vsOthers_longDat[[compareRef]] %>% 
  full_join(ZirFlu_DE_vsOthers_longDat[[compareRef]]) %>%
  filter(valName %in% levels(plotDat_order$valName)) %>%
  mutate(strain = "ZirFlu_H1N1") %>%
  relocate(strain, .before = "compare") %>%
  unite("compare", strain:compare)%>% 
  mutate(compare = gsub("reclassify", "LLvs", compare))


plotDat_order_v2 <- plotDat_order %>%
  full_join(ZirFlu_plotDat) %>%
  mutate(valName = factor(valName, levels = levels(plotDat_order$valName), ordered = TRUE)) %>%
  mutate(compare = factor(compare, 
                          levels = c("B_LLvsLH", "B_LLvsHL", "B_LLvsHH",
                                     "H1N1_LLvsLH", "H1N1_LLvsHL", "H1N1_LLvsHH",
                                     "H3N2_LLvsLH",  "H3N2_LLvsHL",  "H3N2_LLvsHH",
                                     "ZirFlu_H1N1_LLvsLH", "ZirFlu_H1N1_LLvsHL", "ZirFlu_H1N1_LLvsHH")))


plotDat_order_v2 %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order_v3 <- get.plotDat_clusterRow(plotDat_order_v2, colName = "compare", valColumn = "tstatistic")
plotDat_order_v3 %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# LL vs Other (only 1 strain)
DEs_temp <- unique(unlist(res.venn)) # all DEs 
DEs_v2 <- DEs_temp[-which(DEs_temp %in% DEs)] # DEs in only 1 strain

ZirFlu_plotDat <- ZirFlu_tstat_vsOthers_longDat[[compareRef]] %>% 
  full_join(DE_vsOthers_longDat[[compareRef]]) %>%
  filter(valName %in% DEs_v2) %>% filter(strain == "H1N1_reclassify") %>%
  relocate(strain, .before = "compare") %>%
  unite("compare", strain:compare) %>% 
  mutate(compare = gsub("reclassify_reclassify", "LLvs", compare))

plotDat5 <- plotDat_order_iMED %>%
  full_join(ZirFlu_plotDat %>% mutate(compare = paste0("ZirFlu_", compare)))

plotDat5_order <- get.plotDat_clusterRow(plotDat5, 
                                           colName = "compare", valColumn = "tstatistic")
plotDat5_order %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

selected_genes <- c("VEGFA", "TNF", "BTN3A2", "COLEC12", "CLEC7A", "KRT19", "IL6", 
                    "REG4", "PLAUR", "CD160", "PGF", "GZMA", "CSF3", "NCR1", "SPINT2", 
                    "SH2D1A", "CXCL17", "SPINK4", "BCL2L11",
                    "SLAMF7", "ICAM4", "FCAR",
                    "TNFRSF4", "NBN", "CCL20", "CD200R1", "CLEC4G",
                    "FCRL2", "CD22",
                    "CLEC4C", "PARP1", "SIGLEC10",
                    "TGFA", "NPPC", "HLA-E", "SELPLG",
                    "FOXO1", "ITM2A", "NT5C3A", "MAP2K6", "PTX3", "FASLG", "CD79B")
my_custom_lables <- levels(plotDat5_order$valName)
my_custom_lables_v2 <- ifelse(my_custom_lables %in% selected_genes, my_custom_lables, "")


plotDat5_order %>%
  ggplot(aes(x = compare, y = valName, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.y =element_blank()) +
  scale_y_discrete(labels = my_custom_lables_v2)

