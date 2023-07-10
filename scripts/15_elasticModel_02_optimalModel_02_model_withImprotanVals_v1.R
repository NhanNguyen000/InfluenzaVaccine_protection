get_best_result = function(caret_fit) {
  # Aim: get the bet predict models after training models with the data
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

# The rest of the code ====================


## prepare validation-test --------------------------------
valiSets <- list()
# for H1N1
for (year in c("2014", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_T1" = "H1N1_T1")
}

# for H3N2
valiSets$H3N2_2015 <- dat_temp %>% 
  filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_T1" = "H3N2_T1")

# for B
for (year in c("2014", "2015")) {
  valiSets[[paste0("B_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "B_reclassify", "ab_T1" = "B_T1")
}

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_T1" = "Byamagata_T1")


valiSets %>% lapply(function(x) x%>% count(reclassify))


## input variable --------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVariables <- list()
inputVariables$age_sex_proteins <- c("age", "sex", proNames)
inputVariables$age_sex_abT1_proteins <- c("age", "sex", "ab_T1", proNames)
inputVariables$age_sex_abT1_metabolites <- c("age", "sex", "ab_T1", meboNames)
inputVariables$age_sex_abT1_proteins_metabolites <- c("age", "sex", "ab_T1", proNames, meboNames)
inputVariables$age_sex_proteins_metabolites <- c("age", "sex", proNames, meboNames)


inputVal <- "age_sex_abT1_proteins_metabolites"

# prepare input and validation sets
train_inputSet <- trainSet[, inputVariables[[inputVal]]]
train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

validation_inputSets <- valiSets %>% 
  lapply(function(x) x[, inputVariables[[inputVal]]])

validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% 
           as.vector() %>% unlist() %>% as.factor())

## run the model -----------------------------
netFit <- train(x = train_inputSet,
                y = train_predOut,
                method = "glmnet", 
                #metric = 'Accuracy', # for classification
                metric = "Kappa", # similar to classification accuracy but it is useful to normalize the imbalance in classes
                tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                     .lambda = seq(0,1,by=0.01)),
                trControl = trainControl(method="repeatedcv",
                                         number=5,
                                         repeats=10))
netFit_model <- netFit # save model
res_crosVal <- get_best_result(netFit) # Test accuracy

netFit2 <- train(x = train_inputSet,
                 y = train_predOut,
                 method = "glmnet", 
                 #metric = 'Accuracy', # for classification
                 #metric = "Kappa", # similar to classification accuracy but it is useful to normalize the imbalance in classes
                 metric = "ROC",
                 tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                      .lambda = seq(0,1,by=0.01)),
                 trControl = trainControl(method="repeatedcv",
                                          number=5,
                                          repeats=10,
                                          summaryFunction=twoClassSummary, 
                                          classProbs=T,
                                          savePredictions = T),
)

## featurePlot -----------------
# https://topepo.github.io/caret/visualizations.html
featurePlot(x = train_inputSet[, 1:10], 
            y = train_predOut, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))
featurePlot(x = train_inputSet[, 1:4], 
            y = train_predOut, 
            plot = "pairs", 
            auto.key = list(columns = 3))

featurePlot(x = train_inputSet[, 5:9], 
            y = train_predOut, 
            plot = "box",
            scales = list(y = list(relation = "free"),
                          x = list(rot = 90)),
            layout = c(4, 1))

featurePlot(x = train_inputSet[, 1:10], 
            y = as.factor(train_predOut), 
            plot = "box",
            scales = list(y = list(relation = "free")))

# Select important features, using recursive feature elimination (RFE). ----------------
# Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset. 
# Step 2: Keeping priority to the most important variables, iterate through by building models of given subset sizes, that is, subgroups of most important predictors determined from step 1. Ranking of the predictors is recalculated in each iteration. 
# Step 3: The model performances are compared across different subset sizes to arrive at the optimal number and list of final predictors. It can be implemented using the rfe() function and you have the flexibility to control what algorithm rfe uses and how it cross validates by defining the rfeControl().

subsets <- c(1, 5, 10, 15, 20, 30, 40, 50, 100)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x=train_inputSet, 
                 y=as.factor(train_predOut),
                 sizes = subsets,
                 rfeControl = ctrl)

lmProfile
# The top 5 variables (out of 5):
#   ab_T1, SERPINB8, C7H10O6, C20H18O7, C15H24O5
predictors(lmProfile)
plot(lmProfile, type = c("g", "o"))


# See available algorithms in caret
# modelnames <- paste(names(getModelInfo()), collapse=',  ')
# modelnames

netFit_models <- list()
for (i in c(1:5)) {
  netFit_models[[i]] <- train(x = train_inputSet,
                              y = train_predOut,
                              method = "glmnet", 
                              #metric = 'Accuracy', # for classification
                              #metric = "Kappa", # similar to classification accuracy but it is useful to normalize the imbalance in classes
                              metric = "ROC",
                              tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                                   .lambda = seq(0,1,by=0.01)),
                              trControl = trainControl(method="repeatedcv",
                                                       number=5,
                                                       repeats=10,
                                                       summaryFunction=twoClassSummary, 
                                                       classProbs=T,
                                                       savePredictions = T),
  )
}


model_netFit <- netFit_models[[1]]

# elastic net regression parameter tuning
model_netFit$bestTune
ggplot(model_netFit) +
  labs(title = "Elastic Net Regression Parameter Tuning", x = "lambda")

# importantce variables
varimp<- varImp(model_netFit)
plot(varimp, main="Variable Importance")





# check variable distribution --------------------------
library(skimr)
skimmed <- skim_to_wide(train_inputSet)
View(skimmed)

library(pROC)
# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 2
# Plot:
plot.roc(rfFit$pred$obs[selectedIndices],
         rfFit$pred$M[selectedIndices])
