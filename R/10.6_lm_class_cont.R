# 10.6_lm_classification.R

# Purpose:
# Investigating sgRNA classifications made by xgb_model

# Input:

# Import variables and models ---------------------------------------------

source("R/10.5_lm_accuracy.R")
source("R/classification_metrics.R")


# Investigating Training sgRNAs -------------------------------------------

# Levels for observations (negative score = high activity)
training$obs_factor <- ifelse(training$centered_LFC > 0, "low", "high")
# Levels for predictions (negative score = high activity)
training$pred_factor <- ifelse(training$predictions > 0, "low", "high")


# Full Training Confusion Matrix ------------------------------------------

# True Positive, high observed activity predicted as high
train_full_TP <- nrow(training[
  training$obs_factor == "high" & training$pred_factor == "high", ]
)

# True Negative, low observed activity predicted as low
train_full_TN <- nrow(training[
  training$obs_factor == "low" & training$pred_factor == "low",]
  )

# False Postive, low observed activity predicted as high
train_full_FP <- nrow(training[
  training$obs_factor == "low" & training$pred_factor == "high",]
  )
# False Negative, high observed activity predicted as low
train_full_FN <- nrow(training[
  training$obs_factor == "high" & training$pred_factor == "low",]
  )

conf_matrix_tr_full <- matrix(c(train_full_TP, train_full_FP, 
                                train_full_FN, train_full_TN), 
                              nrow=2, byrow=TRUE)
colnames(conf_matrix_tr_full) <- c("Actual high", "Actual low")
rownames(conf_matrix_tr_full) <- c("Predicted high", "Predicted low")
conf_matrix_tr_full

#               Actual high  Actual low
# Predicted high   11684        8612
# Predicted low    8866         12363


# Confusion Matrix on Top/Bottom Training Guides --------------------------

training_sig <- training[training$centered_LFC <= train_q25 | 
                           training$centered_LFC >= train_q75, ]

# True Positive, high observed activity predicted as high
train_sig_TP <- nrow(training_sig[
  training_sig$obs_factor == "high" & training_sig$pred_factor == "high", ]
)

# True Negative, low observed activity predicted as low
train_sig_TN <- nrow(training_sig[
  training_sig$obs_factor == "low" & training_sig$pred_factor == "low",]
)

# False Postive, low observed activity predicted as high
train_sig_FP <- nrow(training_sig[
  training_sig$obs_factor == "low" & training_sig$pred_factor == "high",]
)
# False Negative, high observed activity predicted as low
train_sig_FN <- nrow(training_sig[
  training_sig$obs_factor == "high" & training_sig$pred_factor == "low",]
)

conf_matrix_tr_sig <- matrix(c(train_sig_TP, train_sig_FP, 
                        train_sig_FN, train_sig_TN), 
                      nrow=2, byrow=TRUE)
colnames(conf_matrix_tr_sig) <- c("Actual high", "Actual low")
rownames(conf_matrix_tr_sig) <- c("Predicted high", "Predicted low")
conf_matrix_tr_sig

#               Actual high  Actual low
# Predicted high    6647        3966
# Predicted low     3735        6416        


# Investigating Testing sgRNAs -------------------------------------------

# Levels for observations (negative score = high activity)
testing$obs_factor <- ifelse(testing$centered_LFC > 0, "low", "high")
# Levels for predictions (negative score = high activity)
testing$pred_factor <- ifelse(testing$predictions > 0, "low", "high")


# Full Testing Confusion Matrix ------------------------------------------

# True Positive, high observed activity predicted as high
test_full_TP <- nrow(testing[
  testing$obs_factor == "high" & testing$pred_factor == "high", ]
)

# True Negative, low observed activity predicted as low
test_full_TN <- nrow(testing[
  testing$obs_factor == "low" & testing$pred_factor == "low",]
)

# False Postive, low observed activity predicted as high
test_full_FP <- nrow(testing[
  testing$obs_factor == "low" & testing$pred_factor == "high",]
)
# False Negative, high observed activity predicted as low
test_full_FN <- nrow(testing[
  testing$obs_factor == "high" & testing$pred_factor == "low",]
)

conf_matrix_te_full <- matrix(c(test_full_TP, test_full_FP, 
                        test_full_FN, test_full_TN), 
                      nrow=2, byrow=TRUE)
colnames(conf_matrix_te_full) <- c("Actual high", "Actual low")
rownames(conf_matrix_te_full) <- c("Predicted high", "Predicted low")
conf_matrix_te_full

#              Actual high  Actual low
# Predicted high   2899       2164
# Predicted low    2252       3067
 

# Confusion Matrix on Top/Bottom Training Guides --------------------------

testing_sig <- testing[testing$centered_LFC <= test_q25 | 
                       testing$centered_LFC >= test_q75, ]

# True Positive, high observed activity predicted as high
test_sig_TP <- nrow(testing_sig[
  testing_sig$obs_factor == "high" & testing_sig$pred_factor == "high", ]
)

# True Negative, low observed activity predicted as low
test_sig_TN <- nrow(testing_sig[
  testing_sig$obs_factor == "low" & testing_sig$pred_factor == "low",]
)

# False Postive, low observed activity predicted as high
test_sig_FP <- nrow(testing_sig[
  testing_sig$obs_factor == "low" & testing_sig$pred_factor == "high",]
)
# False Negative, high observed activity predicted as low
test_sig_FN <- nrow(testing_sig[
  testing_sig$obs_factor == "high" & testing_sig$pred_factor == "low",]
)

conf_matrix_te_sig <- matrix(c(test_sig_TP, test_sig_FP, 
                               test_sig_FN, test_sig_TN), 
                      nrow=2, byrow=TRUE)
colnames(conf_matrix_te_sig) <- c("Actual high", "Actual low")
rownames(conf_matrix_te_sig) <- c("Predicted high", "Predicted low")
conf_matrix_te_sig

#               Actual high  Actual low
# Predicted high    1617        999
# Predicted low     979         1597


# Storing Final Results ---------------------------------------------------

get_results_row <- function(conf_matrix, name) {
  df <- data.frame( "Data" = name,
                    "Accuracy" = accuracy(conf_matrix),
                    "Precision" = precision(conf_matrix),
                    "Sensitivity" = sensitivity(conf_matrix),
                    "Specificity" = specificity(conf_matrix),
                    "NPV" = negative_predictive_value(conf_matrix))
  return(df)
}

# Computing evaluation metrics for all confusion matrices
full_tr <- get_results_row(conf_matrix_tr_full, "FullTrainSet")
sig_tr <- get_results_row(conf_matrix_tr_sig, "SigTrainSet")
full_te <- get_results_row(conf_matrix_te_full, "FullTestSet")
sig_te <- get_results_row(conf_matrix_te_sig, "SigTestSet")
results <- rbind(full_tr, sig_tr, full_te, sig_te)

# Adding F1-Score
results$F1Score <- 2*(results$Precision*results$Sensitivity) / 
  (results$Precision + results$Sensitivity)

results

#Data         Accuracy  Precision   Sensitivity Specificity NPV       F1Score
#FullTrainSet 0.5790969 0.5756799   0.5685645   0.5894160   0.5823637 0.5721001
#SigTrainSet  0.6291177 0.6263074   0.6402427   0.6179927   0.6320560 0.6331984
#FullTestSet  0.5746484 0.5725854   0.5628033   0.5863124   0.5766121 0.5676522
#SigTestSet   0.6190293 0.6181193   0.6228814   0.6151772   0.6199534 0.6204912

