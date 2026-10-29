# 14_final_model_exploration.R

# Purpose:
# Explore final linear and xgb models

# Input:
# - Linear model
# - XGB model
# - variables for training, training index and testing index

# Output: 
# - xgb_pred vs lm_pred plots for train and test data -> results/models/

# Reading in necessary libraries ------------------------------------------

library(xgboost)
library(ggplot2)
library(basicstats)

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

final_lm_plot <- ggplot(data = test, 
                       aes(x = observed,
                           y = linear_predictions)) +
  geom_point() +
  labs(title = "Final Linear Model Predictions vs Observations (test data set)",
       x = "Observed Centered Log2 Fold Change",
       y = "Linear Model Prediction") + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw()
ggsave( filename = "results/models/final_lm_fit.png",
        plot = final_lm_plot,
        width = 10, height = 6, dpi = 300)

final_xgb_plot <- ggplot(data = test, 
                        aes(x = observed,
                            y = xgb_predictions)) +
  geom_point() +
  labs(title = "Final XGBoost Model Predictions vs Observations (test data set)",
       x = "Observed Centered Log2 Fold Change",
       y = "XGBoost Model Prediction") + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw()
ggsave( filename = "results/models/final_xgb_fit.png",
        plot = final_xgb_plot,
        width = 10, height = 6, dpi = 300)


# Visualizing Model Similarity --------------------------------------------

# Comparing Models on training data set
xgb_v_lm_tr <- ggplot(data = train, 
                      aes(x=linear_predictions, 
                          y=xgb_predictions, 
                          )) + 
  geom_point(alpha = 0.3) + 
  labs(title = "Train XGBoost vs. Linear",
       x = "Linear",
       y = "XGBoost") + 
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()
ggsave( filename = "results/models/xgb_v_lm_tr.png",
        plot = xgb_v_lm_tr,
        width = 10, height = 6, dpi = 300)

# Comparing Models on test data set
xgb_v_lm_te <- ggplot(data = test, 
                      aes(x=linear_predictions, 
                          y=xgb_predictions, 
                      )) +  
  geom_point(alpha = 0.3) + 
  labs(title = "Test Extreme Gradient Boosting vs. Linear Regression",
       x = "Linear Prediction",
       y = "XGBoost Prediction") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_minimal()
ggsave( filename = "results/models/xgb_v_lm_te.png",
        plot = xgb_v_lm_te,
        width = 10, height = 6, dpi = 300)


# Comparing Predictions ---------------------------------------------------

# Computing residuals on test set predictions
test$lm_res <- abs(test$observed - test$linear_predictions)
test$xgb_res <- abs(test$observed - test$xgb_predictions)

# Computing residuals on training set predictions
train$lm_res <- train$observed - train$linear_predictions
train$xgb_res <- train$observed - train$xgb_predictions


results <- data.frame(corr_test = round(correlation(test$linear_predictions, 
                                              test$xgb_predictions,
                                              "spearman"), 3),
                      corr_train = round(correlation(train$linear_predictions, 
                                               train$xgb_predictions,
                                               "spearman"),3),
                      corr_res_test = round(correlation(test$lm_res,
                                                   test$xgb_res,
                                                   "spearman"), 3),
                      corr_res_train = round(correlation(train$lm_res,
                                                    train$xgb_res,
                                                    "spearman"), 3))
results
#corr_test corr_train corr_res_test corr_res_train
#     0.894      0.891         0.966           0.99

