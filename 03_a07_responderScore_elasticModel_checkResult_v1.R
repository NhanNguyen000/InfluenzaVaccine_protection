load("20230421_res.elasticModel_2groups_1to2Ratio_iMED.RData")
load("20230421_res.elasticModel_2groups_equalRatio_iMED.RData")
load("20230421_res.elasticModel_4groups_iMED.RData")

load("20230421_res.elasticModel_2groups_1to2Ratio_splitCohort.RData")
load("20230421_res.elasticModel_2groups_equalRatio_splitCohort.RData")
# load results from elastic model with using iMED as training cohort and ZirFlu as validation cohort 
res_input <- res.elasticModel_2groups_equalRatio_splitCohort
res_input <- res.elasticModel_2groups_1to2Ratio_splitCohort




# check model results form 2 cohort together------------------------------------------
load("20230421_res.elasticModel_2groups_1to2Ratio.RData")
load("20230421_res.elasticModel_2groups_equalRatio.RData")
load("20230421_res.elasticModel_4groups.RData")

# load results from elastic model with combine cohorts
res_input <- res.elasticModel_4groups
res_input <- res.elasticModel_2groups_equalRatio
res_input <- res.elasticModel_2groups_1to2Ratio


res_crosVal_Tab <- res_input %>%
  lapply(function(x) x$res_crosVal %>% purrr::reduce(rbind))
range(res_crosVal_Tab$proteins$Accuracy)
range(res_crosVal_Tab$age_proteins$Accuracy)
range(res_crosVal_Tab$age_sex_proteins$Accuracy)
range(res_crosVal_Tab$age_sex_abBaseline_proteins$Accuracy)

res_preVal_Tab <- res_input %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(x) x$overall) %>% 
           purrr::reduce(rbind) %>% as.data.frame())
range(res_preVal_Tab$proteins$Accuracy)
range(res_preVal_Tab$age_proteins$Accuracy)
range(res_preVal_Tab$age_sex_proteins$Accuracy)
range(res_preVal_Tab$age_sex_abBaseline_proteins$Accuracy)

res_preVal_Tab_detail <- res_input %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(x) x$byClass) %>% 
           purrr::reduce(rbind) %>% as.data.frame())
modelInfo <- "Sensitivity"
modelType <- "proteins"
modelType <- "age_proteins"
modelType <- "age_sex_proteins"
modelType <- "age_sex_abBaseline_proteins"
#modelInfo <- "Precision"
#modelInfo <- "Recall"

range(res_preVal_Tab_detail[[modelType]][[modelInfo]])

# if res_input is res.elasticModel_4groups
coefs_Tab_all <- res.elasticModel_4groups %>%
  lapply(function(x) x$coefs_Tab %>% 
           lapply(function(x) x %>% column_to_rownames("Protein")) %>%
           imap(~.x %>% rename_with(function(x) paste0(x, "_run_", .y))) %>%
           lapply(function(x) x %>% rownames_to_column("valName")) %>%
           purrr::reduce(full_join) %>% column_to_rownames("valName"))

# if res_input is res.elasticModel for 2 groups
coefs_Tab_all <- res_input %>%
  lapply(function(x) x$coefs_Tab %>% 
           imap(~.x %>% rename_with(function(x) paste0("run_", .y)))%>%
           lapply(function(x) x %>% rownames_to_column("valName")) %>%
           purrr::reduce(full_join) %>% column_to_rownames("valName"))

# check the model performance ----------------------
load("20230404_res.elasticModel_2groups_1to2Ratio_splitCohort_zoomIn.RData")
res.fit <- 
  ggplot(res.fit, aes(x = "")) +
  geom_histogram(binwidth = .05) + facet_wrap(~obs) +
  xlab("Probability of class")

confusionMatrix(data = test_set$pred, reference = test_set$obs)

confusionMatrix(validation_predOut, val.pred) 

res.elasticModel_2groups_1to2Ratio_splitCohort$proteins$res_preVal[[1]]

# the rest to check ---------------
a <- rowMeans(coefs_Tab_all) %>% 
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% 
  rename("coef" = ".")

a2 <- rowMeans(coefs_Tab_all)[which(rowMeans(coefs_Tab_all) != 0)] %>%
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% filter(proteins != "(Intercept)") %>%
  rename("avg_coef" = ".") %>% arrange(avg_coef)

a2 %>%
  ggplot(aes(x = reorder(proteins, avg_coef), y = avg_coef)) + 
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) 

a3 <- a2 %>% mutate("abs_avg_coef" = abs(avg_coef)) %>% arrange(abs_avg_coef)

# only LL 
b <- coefs_Tab_all %>% select(matches("LL"))
ab <- rowMeans(b)  %>% 
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% 
  rename("coef" = ".")

ab2 <- rowMeans(b)[which(rowMeans(b) != 0)] %>%
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% filter("proteins" != "(Intercept)") %>%
  rename("avg_coef" = ".") %>% arrange(avg_coef)

# only LH 
b <- coefs_Tab_all %>% select(matches("LH"))
ab <- rowMeans(b)  %>% 
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% 
  rename("coef" = ".")

ab2 <- rowMeans(b)[which(rowMeans(b) != 0)] %>%
  as.data.frame() %>% 
  rownames_to_column("proteins") %>% filter("proteins" != "(Intercept)") %>%
  rename("avg_coef" = ".") %>% arrange(avg_coef)
