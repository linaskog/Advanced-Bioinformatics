# 10.5_lm_accuracy.R

# Purpose:
# Investigating sgRNA classifications made by xgb_model

# Input:
# - variables
# - train and test indices
# - linear model


# Import variables and models ---------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")

linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds"
)


# Investigating Training sgRNAs -------------------------------------------

# Data frame for storing training metrics
training <- variables[train_index, 
                      c("UID", "Gene", "centered_LFC")]
training$predictions <- predict(linear_model, 
                                variables[train_index, ])

# Comparing observed and predicted scores (positive or negative)
training$same_sign <- sign(training$centered_LFC) == sign(training$predictions)

correct_train_classification <- nrow(training[training$same_sign == TRUE,])
incorrect_train_classification <- nrow(training[training$same_sign == FALSE,])

print(paste0("Number of correct classifications: ",
             correct_train_classification
             )
      )
# Number of correct classifications: 23864

print(paste0("Number of incorrect classifications: ",
             incorrect_train_classification
             )
      )
# Number of incorrect classifications: 17661

print(paste0("Fraction correct classifications: ",
             round(correct_train_classification / nrow(training), 3)
             )
      )
# Fraction correct classifications: 0.575


# Investigating Testing sgRNAs -------------------------------------------

# Data frame for storing testing metrics
testing <- variables[test_index, 
                      c("UID", "Gene", "centered_LFC")]
testing$predictions <- predict(linear_model, 
                                variables[test_index, ])

# Comparing observed and predicted scores (positive or negative)
testing$same_sign <- sign(testing$centered_LFC) == sign(testing$predictions)

correct_test_classification <- nrow(testing[testing$same_sign == TRUE,])
incorrect_test_classification <- nrow(testing[testing$same_sign == FALSE,])

print(paste0("Number of correct classifications: ",
             correct_test_classification
             )
      )
# Number of correct classifications: 5926

print(paste0("Number of incorrect classifications: ",
             incorrect_test_classification
             )
      )
# Number of incorrect classifications: 4456

print(paste0("Fraction correct classifications: ",
             round(correct_test_classification / nrow(testing), 3)
             )
      )
# Fraction correct classifications: 0.571

# Looking Only at Top/Bottom Guides ---------------------------------------

# Training Guides
# Finding top quartile training guide observed scores
train_q75 <- quantile(training$centered_LFC, 0.75)
# Finding bottom quartile training guide observed scores
train_q25 <- quantile(training$centered_LFC, 0.25)

# Data frame with only top and bottom scoring guides
sig_training <- training[training$centered_LFC <= train_q25 | 
                           training$centered_LFC >= train_q75, ]
nrow(sig_training)
# Number of top/bottom ranking genes: 20764

correct_tr_class <- nrow(sig_training[sig_training$same_sign == T, ])
incorrect_tr_class <- nrow(sig_training[sig_training$same_sign == F, ])

print(paste0("Number of correct classifications: ",
             correct_tr_class))
# Number of correct classifications: 13063

print(paste0("Number of incorrect classifications: ",
             incorrect_tr_class))
# Number of incorrect classifications: 7701

print(paste0("Fraction correct classifications: ",
             round(correct_tr_class / nrow(sig_training), 3)))
# Fraction correct classifications: 0.629

# Testing Guides
# Finding top quartile testing guide observed scores
test_q75 <- quantile(testing$centered_LFC, 0.75)
# Finding bottom quartile training guide observed scores
test_q25 <- quantile(testing$centered_LFC, 0.25)

# Data frame with only top and bottom scoring guides
sig_testing <- testing[testing$centered_LFC <= test_q25 | 
                           testing$centered_LFC >= test_q75, ]

correct_te_class <- nrow(sig_testing[sig_testing$same_sign == T, ])
incorrect_te_class <- nrow(sig_testing[sig_testing$same_sign == F, ])

print(paste0("Number of correct classifications: ",
             correct_te_class))
# Number of correct classifications: 3214

print(paste0("Number of incorrect classifications: ",
             incorrect_te_class))
# Number of incorrect classifications: 1978

print(paste0("Fraction correct classifications: ",
             round(correct_te_class / nrow(sig_testing), 3)))
# Fraction correct classifications: 0.619
