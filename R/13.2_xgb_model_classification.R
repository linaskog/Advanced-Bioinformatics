# 13.2_xgb_model_classification.R

# Purpose:
# Investigating sgRNA classifications made by xgb_model

# Input:
# Sources R/13.1_xgb_model_genes.R


# Import variables and models ---------------------------------------------

source("R/13.1_xgb_model_genes.R")

# Investigating Training sgRNAs -------------------------------------------

# Comparing observed and predicted scores (positive or negative)
training$same_sign <- sign(training$centered_LFC) == sign(training$predictions)

# Correct classifications if same sign for observed and predicted scores
correct_train_classification <- nrow(training[training$same_sign == TRUE,])
incorrect_train_classification <- nrow(training[training$same_sign == FALSE,])

print(paste0("Number of correct classifications: ",
             correct_train_classification
             )
      )
# Number of correct classifications: 23984

print(paste0("Number of incorrect classifications: ",
             incorrect_train_classification
             )
      )
# Number of incorrect classifications: 17541

print(paste0("Fraction correct classifications: ",
             round(correct_train_classification / nrow(training), 3)
             )
      )
# Fraction correct classifications: 0.578


# Investigating Testing sgRNAs -------------------------------------------

# Comparing observed and predicted scores (positive or negative)
testing$same_sign <- sign(testing$centered_LFC) == sign(testing$predictions)

# Correct classifications if same sign for observed and predicted scores
correct_test_classification <- nrow(testing[testing$same_sign == TRUE,])
incorrect_test_classification <- nrow(testing[testing$same_sign == FALSE,])

print(paste0("Number of correct classifications: ",
             correct_test_classification)
      )
# Number of correct classifications: 5958

print(paste0("Number of incorrect classifications: ",
             incorrect_test_classification)
      )
# Number of incorrect classifications: 4424

print(paste0("Fraction correct classifications: ",
             round(correct_test_classification / nrow(testing), 3))
      )
# Fraction correct classifications: 0.574


# Looking Only at Top/Bottom Guides ---------------------------------------

# Training Set ------------------------------------------------------------

# Finding top quartile training guide observed scores
train_q75 <- quantile(training$centered_LFC, 0.75)
# Finding bottom quartile training guide observed scores
train_q25 <- quantile(training$centered_LFC, 0.25)

# Data frame with only top and bottom scoring guides
sig_training <- training[training$centered_LFC <= train_q25 | 
                           training$centered_LFC >= train_q75, ]
nrow(sig_training)
# Number of significant guides (0.25, 0.75): 20764

correct_tr_class <- nrow(sig_training[sig_training$same_sign == T, ])
incorrect_tr_class <- nrow(sig_training[sig_training$same_sign == F, ])

print(paste0("Number of correct classifications: ",
             correct_tr_class))
# Number of correct classifications: 13060

print(paste0("Number of incorrect classifications: ",
             incorrect_tr_class))
# Number of incorrect classifications: 7704

print(paste0("Fraction correct classifications: ",
             round(correct_tr_class / nrow(sig_training), 3)))
# Fraction correct classifications: 0.629


# Testing Set -------------------------------------------------------------

# Finding top quartile testing guide observed scores
test_q75 <- quantile(testing$centered_LFC, 0.75)
# Finding bottom quartile training guide observed scores
test_q25 <- quantile(testing$centered_LFC, 0.25)

# Data frame with only top and bottom scoring guides
sig_testing <- testing[testing$centered_LFC <= test_q25 | 
                         testing$centered_LFC >= test_q75, ]
nrow(sig_testing)
# Number of significant guides (0.25, 0.75): 5192

correct_te_class <- nrow(sig_testing[sig_testing$same_sign == T, ])
incorrect_te_class <- nrow(sig_testing[sig_testing$same_sign == F, ])

print(paste0("Number of correct classifications: ",
             correct_te_class))
# Number of correct classifications: 3231

print(paste0("Number of incorrect classifications: ",
             incorrect_te_class))
# Number of incorrect classifications: 1961

print(paste0("Fraction correct classifications: ",
             round(correct_te_class / nrow(sig_testing), 3)))
# Fraction correct classifications: 0.622
