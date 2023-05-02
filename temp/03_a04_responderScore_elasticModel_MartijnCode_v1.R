rm(list = ls())

library(dplyr)
library(magrittr)
library(caret)
library(openxlsx)

setwd('/vol/projects/mzoodsma/prediction_immunological_age/')
load('data/500FG.RData')
load('data/300BCG.RData')

# Only common celltypes
conv <- read.xlsx('celltypes_conversion.xlsx') %>% filter(include)
FG500 <- FG500[, conv$FG]
BCG300 <- BCG300[, conv$BCG_match]

colnames(FG500) <- colnames(BCG300)

prep_model <- preProcess(FG500)

train_x <- predict(prep_model, FG500)
train_y <- meta.FG500

validation_x <- predict(prep_model, BCG300)
validation_y <- meta.bcg300$age

n <- 50

test.res <- list()
val.res <- list()
coefficients <- list()

test.preds <- list()
val.preds <- list()

for (i in 1:n) {
  netFit <- train(x = train_x,
                  y = train_y,
                  method = "glmnet", metric = 'RMSE',
                  tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                       .lambda = seq(0,1,by=0.01)),
                  trControl = trainControl(method="repeatedcv",
                                           number=5,
                                           repeats=10))
  
  # Test accuracy
  test.pred <- predict(netFit, newdata = train_x, type = 'raw')
  res1 <- cor.test(train_y, test.pred)
  test.res[[i]] <- data.frame(
    'cor' = res1$estimate, 
    'pval' = res1$p.value
  )
  
  # Training accuracy
  val.pred <- predict(netFit, newdata = validation_x, type = 'raw')
  res2 <- cor.test(validation_y, val.pred)
  val.res[[i]] <- data.frame(
    'cor' = res2$estimate, 
    'pval' = res2$p.value
  )
  
  
  # Coefficients
  coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
  coefs <- data.frame(coefs[, 1])
  colnames(coefs) <- paste0('run_', i)
  coefficients[[i]] <- coefs
  
  # Predictions on the test set and the validation set. Mean and SD
  test.preds[[i]] <- test.pred
  val.preds[[i]] <- val.pred
  
}
