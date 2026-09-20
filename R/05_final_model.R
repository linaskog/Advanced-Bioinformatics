#05_FINAL MODEL.R

# Reading in Libraries ----------------------------------------------------

library(xgboost)

# Reading in Data ---------------------------------------------------------

source("R/04_model_optimization.R")


# Functions ---------------------------------------------------------------

correlation <- function(obs, pred, method){
  return( cor(pred, obs, method = method))
}

RMSE <- function(obs, pred){
  return( sqrt(mean((obs-pred)^2)))
}

MAE <- function(obs, pred){
  return(mean(abs(obs-pred)))
}

compile <- function(model, rounds) {
  pred_train <- predict(model, dtrain)
  pred_test <- predict(model, dtest)
  
  data.frame(
    rounds = rounds,
    trainRMSE = RMSE(y[train_index], pred_train),
    testRMSE = RMSE(y[test_index], pred_test),
    trainMAE = MAE(y[train_index], pred_train),
    testMAE = MAE(y[test_index], pred_test),
    trainSpearman = correlation(y[train_index], pred_train, "spearman"),
    testSpearman = correlation(y[test_index], pred_test, "spearman")
  )
}


# Training Final Model ----------------------------------------------------

xgb_model <- xgb.train(
  params = final_params, 
  data = dtrain,
  nrounds = final_nrounds)

# Final Model Metrics

xgb_model

attributes(xgb_model)

final_results <- compile(xgb_model, final_cv$rounds)

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
