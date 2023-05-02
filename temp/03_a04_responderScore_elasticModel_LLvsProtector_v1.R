library(dplyr)
library(magrittr)
library(caret)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

# input data --------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# impute protein
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay") %>% drop_na()

# use 2 group LL vs. Protector (accuracy is higher 0.8, but it is because everything is assign for protector) ----------
dat <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))

dat %>% count(reclassify) # LL group = 48 people
dat_v2 <- dat %>% group_by(reclassify) %>% sample_n(size = 48)
dat_v2 %>% count(reclassify) #equal LL (n = 48) and protector (n = 48)

## set the seed to make your partition reproducible
set.seed(123)
trainSet <- dat_v2 %>% group_by(reclassify) %>% slice_sample(prop = 0.70)
trainSet %>% count(reclassify) # per group n = 33

validateSet <- dat_v2 %>% anti_join(trainSet)
validateSet %>% count(reclassify) # per group n = 15

# run the model
train_x1 <- trainSet[, rownames(inputDat)]
train_y1 <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

validation_x1 <- validateSet[, rownames(inputDat)]
validation_y1 <- validateSet[, c("reclassify")] %>% as.vector() %>% unlist() %>% as.factor()

n <- 50
res_crosVal <- list()
res_preVal <- list()
coefs_Tab <- list()
for (i in 1:n) {
  netFit <- train(x = train_x1,
                  y = train_y1,
                  method = "glmnet", metric = 'Accuracy', # for classification
                  tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                       .lambda = seq(0,1,by=0.01)),
                  trControl = trainControl(method="repeatedcv",
                                           number=5,
                                           repeats=10))
  res_crosVal[[i]] <- get_best_result(netFit) # Test accuracy
  
  
  val.pred <- predict(netFit, newdata = validation_x1, type = 'raw')
  res_preVal[[i]] <- confusionMatrix(validation_y1, val.pred) # Training accuracy
  
  # Coefficients
  coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
  coefs_Tab[[i]] <- coefs[, 1] %>% as.data.frame()
  print(paste0("Run ", i, " time."))
}

# save(res_crosVal, res_preVal, coefs_Tab , file = "elasticModel_LLvsProtector.RData")
load("elasticModel_LLvsProtector.RData")
# check results
res_crosVal_Tab <- res_crosVal %>% purrr::reduce(rbind)
range(res_crosVal_Tab$Accuracy)

res_preVal_Tab <- res_preVal %>% lapply(function(x) x$overall) %>% 
  purrr::reduce(rbind) %>% as.data.frame()
range(res_preVal_Tab$Accuracy)

coefs_Tab_all <- coefs_Tab %>% 
  imap(~.x %>% rename_with(function(x) paste0("run_", .y))) %>%
  lapply(function(x) x %>% rownames_to_column("Proteins")) %>%
  purrr::reduce(full_join) %>% column_to_rownames("Proteins")

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
