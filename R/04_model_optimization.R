# 04_model_optimization.R
# Advanced Bioinformatics: sgRNA quality project

# Purpose:
# Model sgRNA efficacy on sequence composition

# Input:
# DMatrices containing training and testing data from 03_data_preprocessing.R

# Reading in Libraries ----------------------------------------------------

library(xgboost)


# Reading in Data ---------------------------------------------------------

source("R/03_data_preprocessing.R")

# XGBoost -----------------------------------------------------------------
# nrounds:
# params:
### objective
### nthreads
### max_depth
### eta/learning rate
###

model <- xgb.train(
  data = dtrain,
  nrounds = 1000,
  params = xgb.params(
    eta = 0.04,
    gamma = 0,
    max_depth = 6,
    min_child_weight = 10,
    max_delta_step = 0,
    sampling_method = "uniform",
    objective = "reg:squarederror",
    subsample = 0.8,
    colsample_bytree = 0.5,
    colsample_bylevel = 0.5,
    colsample_bynode = 0.5,
    lambda = 1,
    alpha = 0,
    tree_method = "auto"
  )
)

# Predicting with model

pred_train <- predict(model, dtrain)
pred_test <- predict(model, dtest)

plot(y[train_index], pred_train, main = "Predicted vs. Observed",
     xlab = "Observed", ylab = "Predicted") + abline(1,1, col = "red")

plot(y[test_index], pred_test, main = "Predicted vs. Observed",
     xlab = "Observed", ylab = "Predicted") + abline(1,1, col = "red")

cor(y[train_index], pred_train, method = "spearman")
cor(y[test_index], pred_test, method = "spearman")

#xgb.importance(model)
