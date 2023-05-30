rm(list = ls())

library(tidyverse)
library(caret)
# load the data -------------------------
load("20230522_res.elasticModel.RData")

# load results from elastic model with combine cohorts



res_crosVal_Tab <- res.elasticModel %>%
  lapply(function(x) x$res_crosVal %>% purrr::reduce(rbind))
res_crosVal_Tab %>% lapply(function(x) range(x$Accuracy))

res_preVal_Tab <- res.elasticModel %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(y) y %>% lapply(function(z) z$overall) %>% 
                    purrr::reduce(rbind) %>% as.data.frame()))
res_preVal_Tab %>% 
  lapply(function(x) x %>% lapply(function(y) range(y$Accuracy)))
res_preVal_Tab$age_sex_proteins %>% lapply(function(x) range(x$Accuracy))
res_preVal_Tab$age_sex_abT1_proteins_metabolites %>% lapply(function(x) range(x$Accuracy))
res_preVal_Tab$age_sex_proteins_metabolites %>% lapply(function(x) range(x$Accuracy))

res_preVal_Tab_v2 <- res.elasticModel %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(y) y %>% lapply(function(z) z$byClass) %>% 
                    purrr::reduce(rbind) %>% as.data.frame()))
res_preVal_Tab_v2$age_sex_abT1_proteins_metabolites %>% lapply(function(x) range(x$Recall))
res_preVal_Tab_v2$age_sex_proteins_metabolites %>% lapply(function(x) range(x$Recall))

# compare model performances -------------
https://machinelearningmastery.com/compare-models-and-select-the-best-using-the-caret-r-package/
https://easystats.github.io/performance/

https://bookdown.org/gmli64/do_a_data_science_project_in_10_days/multiple-models-comparison.html
https://machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/
https://towardsdatascience.com/a-guide-to-using-caret-in-r-71dec0bda208

# check variable improtant
https://towardsdatascience.com/create-predictive-models-in-r-with-caret-12baf9941236

# show some plot ---------
twoClassSummary()

set.seed(144)
true_class <- factor(sample(paste0("Class", 1:2), 
                            size = 1000,
                            prob = c(.2, .8), replace = TRUE))
true_class <- sort(true_class)
class1_probs <- rbeta(sum(true_class == "Class1"), 4, 1)
class2_probs <- rbeta(sum(true_class == "Class2"), 1, 2.5)
test_set <- data.frame(obs = true_class,
                       Class1 = c(class1_probs, class2_probs))
test_set$Class2 <- 1 - test_set$Class1
test_set$pred <- factor(ifelse(test_set$Class1 >= .5, "Class1", "Class2"))
ggplot(test_set, aes(x = Class1)) + 
  geom_histogram(binwidth = .05) + 
  facet_wrap(~obs) + 
  xlab("Probability of Class #1")

confusionMatrix(data = test_set$pred, reference = test_set$obs)
View(res.elasticModel$age_sex_proteins$netFit_model)

twoClassSummary(test_set, lev = levels(test_set$obs))
test_set2 <- res_input$age_sex_proteins$res_preVal$H1N1_2014[[1]]