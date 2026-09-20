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

tuning <- function(seed = 2026,
                   data = dtrain,
                   eta = 0.3, # learning rate (0-1)
                   gamma = 0, # minimum loss reduction to make split (0-inf)
                   max_depth = 6, # maximum depth of each decision tree, complexity of model
                   subsample = 1, # proportion training rows randomly sampled to grow each tree
                   colsample_bytree = 1, # ratio of features randomly selected to build each tree
                   min_child_weight = 1
){
  set.seed(seed)
  param_list = list(objective = "reg:squarederror", 
                    eta = eta, 
                    gamma = gamma, 
                    max_depth = max_depth, 
                    subsample = subsample, 
                    colsample_bytree = colsample_bytree, 
                    min_child_weight = min_child_weight)
  
  xgbcv = xgb.cv(params = param_list,
                 data = data,
                 nrounds = 1000,
                 nfold = 5,
                 print_every_n = 10,
                 early_stopping_rounds = 30,
                 maximize = F,
                 metrics = c("rmse", "mae"))
  
  no_rounds <- which.min(xgbcv$evaluation_log$test_rmse_mean)
  
  return(
    list(rounds = no_rounds,
         cvRMSE = xgbcv$evaluation_log$test_rmse_mean[no_rounds],
         cvMAE = xgbcv$evaluation_log$test_mae_mean[no_rounds]))
}

# Tuning Learning Rate ----------------------------------------------------


eta_tune = c(0.01, 0.05, 0.1, 0.15, 0.3, 0.4, 0.5, 0.6)

eta_results <- data.frame(eta = eta_tune,
                          rounds = NA,
                          cvRMSE = NA,
                          cvMAE = NA)

for (i in 1:nrow(eta_results)) {
  x <- tuning(eta = eta_results$eta[i])
  eta_results$rounds[i] <- x$rounds
  eta_results$cvRMSE[i] <- x$cvRMSE
  eta_results$cvMAE[i] <- x$cvMAE
}

best_eta <- eta_results$eta[which.min(eta_results$cvRMSE)]


# Tuning Gamma ------------------------------------------------------------

gamma_tune = c(0, 0.05, 0.1, 0.15, 0.2, 0.3)

gamma_results <- data.frame(gamma = gamma_tune,
                          rounds = NA,
                          cvRMSE = NA,
                          cvMAE = NA)

for (i in 1:nrow(gamma_results)) {
  x <- tuning(eta = best_eta, gamma = gamma_results$gamma[i])
  gamma_results$rounds[i] <- x$rounds
  gamma_results$cvRMSE[i] <- x$cvRMSE
  gamma_results$cvMAE[i] <- x$cvMAE
}

best_gamma <- gamma_results$gamma[which.min(gamma_results$cvRMSE)]

# Tuning max_depth ----------------------------------------------------------

max_depth_tune = c(3, 4, 5, 6, 7, 8, 9)

max_depth_results <- data.frame(max_depth = max_depth_tune,
                          rounds = NA,
                          cvRMSE = NA,
                          cvMAE = NA)

for (i in 1:nrow(max_depth_results)) {
  x <- tuning(eta = best_eta, gamma = best_gamma, max_depth = max_depth_results$max_depth[i])
  max_depth_results$rounds[i] <- x$rounds
  max_depth_results$cvRMSE[i] <- x$cvRMSE
  max_depth_results$cvMAE[i] <- x$cvMAE
}

best_max_depth <- max_depth_results$max_depth[which.min(max_depth_results$cvRMSE)]


# Tuning Subsample --------------------------------------------------------

subsample_tune = c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

subsample_results <- data.frame(subsample = subsample_tune,
                            rounds = NA,
                            cvRMSE = NA,
                            cvMAE = NA)

for (i in 1:nrow(subsample_results)) {
  x <- tuning(eta = best_eta, 
              gamma = best_gamma, 
              max_depth = best_max_depth,
              subsample = subsample_results$subsample[i])
  subsample_results$rounds[i] <- x$rounds
  subsample_results$cvRMSE[i] <- x$cvRMSE
  subsample_results$cvMAE[i] <- x$cvMAE
}

best_subsample <- subsample_results$subsample[which.min(subsample_results$cvRMSE)]


# Tuning colsample_bytree -------------------------------------------------

colsample_bytree_tune = c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

colsample_bytree_results <- data.frame(colsample_bytree = colsample_bytree_tune,
                                rounds = NA,
                                cvRMSE = NA,
                                cvMAE = NA)

for (i in 1:nrow(colsample_bytree_results)) {
  x <- tuning(eta = best_eta, 
              gamma = best_gamma, 
              max_depth = best_max_depth,
              subsample = best_subsample,
              colsample_bytree = colsample_bytree_results$colsample_bytree[i])
  colsample_bytree_results$rounds[i] <- x$rounds
  colsample_bytree_results$cvRMSE[i] <- x$cvRMSE
  colsample_bytree_results$cvMAE[i] <- x$cvMAE
}

best_colsample_bytree <- colsample_bytree_results$colsample_bytree[
  which.min(colsample_bytree_results$cvRMSE)]



# Tuning min_child_weight -------------------------------------------------

min_child_weight_tune = c(0, 0.5, 1, 2, 5, 10)

min_child_weight_results <- data.frame(min_child_weight = min_child_weight_tune,
                                       rounds = NA,
                                       cvRMSE = NA,
                                       cvMAE = NA)

for (i in 1:nrow(min_child_weight_results)) {
  x <- tuning(eta = best_eta, 
              gamma = best_gamma, 
              max_depth = best_max_depth,
              subsample = best_subsample,
              colsample_bytree = best_colsample_bytree,
              min_child_weight = min_child_weight_results$min_child_weight[i])
  min_child_weight_results$rounds[i] <- x$rounds
  min_child_weight_results$cvRMSE[i] <- x$cvRMSE
  min_child_weight_results$cvMAE[i] <- x$cvMAE
}

best_min_child_weight <- min_child_weight_results$min_child_weight[
  which.min(min_child_weight_results$cvRMSE)]


# Final Model Parameters ------------------------------------------------------

final_cv <- tuning(
  eta = best_eta,
  gamma = best_gamma, 
  max_depth = best_max_depth, 
  subsample = best_subsample, 
  colsample_bytree = best_colsample_bytree, 
  min_child_weight = best_min_child_weight)

final_params <- list(objective = "reg:squarederror", 
                     eta = best_eta, 
                     gamma = best_gamma, 
                     max_depth = best_max_depth, 
                     subsample = best_subsample, 
                     colsample_bytree = best_colsample_bytree, 
                     min_child_weight = best_min_child_weight)

final_nrounds <- final_cv$rounds












