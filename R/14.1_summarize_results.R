# 14.1_summarize_results.R

# Purpose:
# Compute summary tables for summarizing findings

# Input:
# - Sources "R/14_final_model_exploration.R",
# - sources "R/13.3_model_class_cont.R"
# - sources "R/10.6_lm_class_cont.R"


# Sourcing Prediction -----------------------------------------------------

source("R/14_final_model_exploration.R")


# Regression Metrics Summary ----------------------------------------------

summary <- data.frame(
  "Model" = c("Linear Model", "XGBoost Model"),
  "TrainR2" = c(round(R2(train$observed, train$linear_predictions),3),
                round(R2(train$observed, train$xgb_predictions),3)),
  "TestR2" = c(round(R2(test$observed, test$linear_predictions),3),
               round(R2(test$observed, test$xgb_predictions),3)),
  "TrainRMSE" = c(round(RMSE(train$observed, train$linear_predictions), 3),
                  round(RMSE(train$observed, train$xgb_predictions), 3)),
  "TestRMSE" = c(round(RMSE(test$observed, test$linear_predictions),3),
                 round(RMSE(test$observed, test$xgb_predictions), 3))
)

summary

#          Model TrainR2 TestR2 TrainRMSE TestRMSE
# Linear Model   0.058  0.046     0.323    0.327
# XGBoost Model   0.064  0.050     0.322    0.326

# Top Coefficient Summary -------------------------------------------------

# Linear Model Coefficients
lm_coeff <- summary(linear_model)$coefficients
lm_coeff <- lm_coeff[order(lm_coeff[, "Pr(>|t|)"]),]
lm_coeff <- lm_coeff[rownames(lm_coeff) != "(Intercept)", ]
top_lm_coeff <- rownames(lm_coeff)[1:5]

# XGBoost Model
shap_values <- shap.values(
  xgb_model = xgb_model,
  X_train = X[train_index, ]
)
shap_scores <- as.matrix(shap_values$mean_shap_score)
shap_scores <- as.matrix(shap_scores[order(-shap_scores), ])
top_xgb_coeff <- rownames(shap_scores)[1:5]

coeff_summary <- data.frame(
  "Linear" <- top_lm_coeff,
  "XGBoost" <- top_xgb_coeff
)

coeff_summary
# Linear              XGBoost
# pos19_G             pos20_G
# pos20_G             pos20_C
# pos19_G:pos20_C     pos17_T
# pos19_C             pos18_C
# pos18_C             pos19_T


# Classification Metrics Summary ------------------------------------------

# XGBoost Model
source("R/13.3_xgb_model_class_cont.R")
# Loading confusion matrix for xgboost model, training data all guides
conf_matrix_tr_xgb <- conf_matrix_tr_full
conf_matrix_tr_xgb
#               Actual high Actual low
#Predicted high       11784       8585
#Predicted low         8766      12390

# Loading confusion matrix for xgboost model, testing set all guides
conf_matrix_te_xgb <- conf_matrix_te_full
conf_matrix_te_xgb
#                   Actual high Actual low
# Predicted high        2916       2150
# Predicted low         2235       3081

# Linear Model
source("R/10.6_lm_class_cont.R")
# Loading confusion matrix for linear model, training data all guides
conf_matrix_tr_lm <- conf_matrix_tr_full
conf_matrix_tr_lm
#                Actual high Actual low
# Predicted high       11684       8612
# Predicted low         8866      12363

conf_matrix_te_lm <- conf_matrix_te_full
conf_matrix_te_lm
#                   Actual high Actual low
# Predicted high        2899       2164
# Predicted low         2252       3067

# Function for computing classification metrics
get_results_row <- function(conf_matrix, name) {
  df <- data.frame( "Data" = name,
                    "Accuracy" = round(basicstats::accuracy(conf_matrix), 3),
                    "Precision" = round(precision(conf_matrix), 3),
                    "Sensitivity" = round(sensitivity(conf_matrix), 3),
                    "Specificity" = round(specificity(conf_matrix), 3),
                    "NPV" = round(NPV(conf_matrix), 3))
  
  return(df)
}

lm_tr_metrics <- get_results_row(conf_matrix_tr_lm, "LMTrain")
lm_te_metrics <- get_results_row(conf_matrix_te_lm, "LMTest")
xgb_tr_metrics <- get_results_row(conf_matrix_tr_xgb, "XGBoostTrain")
xgb_te_metrics <- get_results_row(conf_matrix_te_xgb, "XGBoostTest")
class_metrics <- rbind(lm_tr_metrics,
                       lm_te_metrics,
                       xgb_tr_metrics,
                       xgb_te_metrics)

class_metrics$F1Score <- round(2*(class_metrics$Precision*class_metrics$Sensitivity) / 
                                 (class_metrics$Precision + class_metrics$Sensitivity), 3)
class_metrics

#         Data Accuracy Precision Sensitivity Specificity   NPV F1Score
#      LMTrain    0.579     0.576       0.569       0.589 0.582   0.572
#       LMTest    0.575     0.573       0.563       0.586 0.577   0.568
# XGBoostTrain    0.582     0.579       0.573       0.591 0.586   0.576
#  XGBoostTest    0.578     0.576       0.566       0.589 0.580   0.571

