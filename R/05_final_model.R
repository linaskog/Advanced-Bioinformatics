#05_FINAL MODEL.R

# Reading in Libraries ----------------------------------------------------

library(xgboost)
library(SHAPforxgboost)
#install.packages("basicstats", repos = NULL, type = "source")
library(basicstats)

# Reading in Data ---------------------------------------------------------

source("R/03_data_preprocessing.R")

final_params <- readRDS("final_params.rds")

final_param_list <- list(
  objective = "reg:squarederror",
  eta = final_params$eta,
  gamma = final_params$gamma,
  max_depth = final_params$max_depth,
  subsample = final_params$subsample,
  colsample_bytree = final_params$colsample_bytree,
  min_child_weight = final_params$min_child_weight
)


# Functions ---------------------------------------------------------------


compile <- function(model, rounds) {
  pred_train <- predict(model, dtrain)
  pred_test <- predict(model, dtest)
  
  data.frame(
    rounds = rounds,
    trainRMSE = RMSE(y[train_index], pred_train),
    testRMSE = RMSE(y[test_index], pred_test),
    trainMAE = MAE(y[train_index], pred_train),
    testMAE = MAE(y[test_index], pred_test),
    trainR2 = R2(y[train_index], pred_train),
    testR2 = R2(y[test_index], pred_test),
    trainSpearman = correlation(y[train_index], pred_train, "spearman"),
    testSpearman = correlation(y[test_index], pred_test, "spearman")
  )
}


# Training Final Model ----------------------------------------------------

xgb_model <- xgb.train(
  params = final_param_list, 
  data = dtrain,
  nrounds = final_params$best_rounds)

# Final Model Metrics

xgb_model

attributes(xgb_model)

final_results <- compile(xgb_model, final_params$best_rounds)

# Plotting Variable Importance

var_imp <- xgb.importance(
  feature_names = colnames(dtrain),
  model = xgb_model
)
xgb.plot.importance(var_imp)


# Predicted vs. Observed Plots --------------------------------------------

pred_train <- predict(xgb_model, dtrain)
pred_test <- predict(xgb_model, dtest)

plot(y[train_index], pred_train , main = "Training Set Predicted vs. Observed",
     xlab = "Observed", ylab = "Predicted")
abline(0,1, col = "red")

plot(y[test_index], pred_test, main = "Testing Set Predicted vs. Observed",
     xlab = "Observed", ylab = "Predicted")
abline(0,1, col = "red")


# SHAP analysis ----------------------------------------------------------

shap_values <- shap.values(
  xgb_model = xgb_model,
  X_train = X[train_index,]
)

shap_values$mean_shap_score


shap_long <- shap.prep(shap_contrib = shap_values$shap_score, 
                       X_train = X[train_index,])

shap.plot.summary(shap_long, dilute = 100)
