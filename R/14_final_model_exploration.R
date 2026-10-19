# 14_final_model_exploration.R

# Purpose:
# Explore final linear and xgb models

# Input:
# - Linear model
# - XGB model
# - variables for training, training index and testing index


# Reading in necessary libraries ------------------------------------------

library(xgboost)
library(ggplot2)

# Reading in Models -------------------------------------------------------

xgb_model <- xgb.load(
  "~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")
linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")


# Reading in Data ---------------------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")
train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")



# Preparing data ----------------------------------------------------------

pos <- as.matrix(variables[, grep("pos", colnames(variables))])
additional_data <- as.matrix(variables[, c("GC", "strand", "chromosome", "region_type",
                                           "gene_overlap_count", "H3K27ac_peak_overlap")
])

X <- cbind(pos, additional_data)

# Data frame for testing data
test <- data.frame(
  sgRNA_id = variables[test_index, "UID"],
  guide = variables[test_index, "Guide"],
  observed = variables[test_index, "centered_LFC"],
  linear_predictions = predict(linear_model, variables[test_index,]),
  xgb_predictions = predict(xgb_model, X[test_index,])
)

# Data frame for training data
train <- data.frame(sgRNA_id = variables[train_index, "UID"],
                    guide = variables[train_index, "Guide"],
                    observed = variables[train_index, "centered_LFC"],
                    linear_predictions = predict(linear_model,
                                                 variables[train_index,]),
                    xgb_predictions = predict(xgb_model, 
                                              X[train_index,]))


# Visualizing Model Similarity --------------------------------------------

# Comparing Models on training data set
ggplot(data = train, aes(x=linear_predictions, 
                         y=xgb_predictions, 
                         alpha = 0.3)) + 
  geom_point() + 
  labs(title = "Train XGBoost vs. Linear",
       x = "Linear",
       y = "XGBoost") + 
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()

# Comparing Models on test data set
ggplot(data = test, aes(x=linear_predictions, 
                        y=xgb_predictions,
                        alpha = 0.3)) + 
  geom_point() + 
  labs(title = "Test XGBoost vs. Linear",
       x = "Linear",
       y = "XGBoost") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()


# Comparing Predictions ---------------------------------------------------

# Computing residuals on test set predictions
test$lm_res <- abs(test$observed - test$linear_predictions)
test$xgb_res <- abs(test$observed - test$xgb_predictions)

# Computing residuals on training set predictions
train$lm_res <- train$observed - train$linear_predictions
train$xgb_res <- train$observed - train$xgb_predictions

results <- data.frame(corr_test = correlation(test$linear_predictions, 
                                              test$xgb_predictions),
                      corr_train = correlation(train$linear_predictions, 
                                               train$xgb_predictions),
                      corr_res_test = correlation(test$lm_res,
                                                   test$xgb_res),
                      corr_res_train = correlation(train$lm_res,
                                                    train$xgb_res))
results
#   corr_test corr_train corr_res_test corr_res_train
#   0.8926498  0.8920192     0.9872019      0.9935959
