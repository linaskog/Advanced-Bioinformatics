# 10_linear_model_training.R
# Advanced Bioinformatics: sgRNA quality project

# Purpose:
# Model sgRNA efficacy using sequence composition, reference genome annotation, 
# and ChIP-Seq annotation data

# Input: Variable data frame (R/variables.rds from 09_data_preparation.R)
#
# Output: Linear Model saved as RDS object
#         test and training indexes for model training as RDS


# Loading Libraries -------------------------------------------------------

library(basicstats)

# Load Necessary Data -----------------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

# Splitting Data ----------------------------------------------------------

set.seed(2026)

# Splitting dataset into testing and training sgRNAs based on gene names
genes <- unique(variables$Gene)

train_genes <- sample(genes, floor(0.8 * length(genes)))
test_genes <- setdiff(genes, train_genes)

train_index <- which(variables$Gene %in% train_genes)
test_index <- which(variables$Gene %in% test_genes)

# Saving Train and Test Indeces For Downstream Modeling
saveRDS(train_index, 
        "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
saveRDS(test_index, 
        "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")


# Function for Summarizing Results ----------------------------------------

get_summary <- function(model1, model2, name1, name2) {
  
  y <- variables[, "centered_LFC"]
  
  train_observations <- y[train_index]
  test_observations <- y[test_index]
  
  model1_train_pred <- predict(model1, variables[train_index, ])
  model1_test_pred <- predict(model1, variables[test_index, ])
  model2_train_pred <- predict(model2, variables[train_index, ])
  model2_test_pred <- predict(model2, variables[test_index, ])
  
  return(data.frame( "model" = c(name1, name2),
                     "trainRMSE" = c(RMSE(train_observations, model1_train_pred),
                                     RMSE(train_observations, model2_train_pred)),
                     "testRMSE" = c(RMSE(test_observations, model1_test_pred),
                                    RMSE(test_observations, model2_test_pred)),
                     "trainMAE" = c(MAE(train_observations, model1_train_pred),
                                    MAE(train_observations, model2_train_pred)),
                     "testMAE" = c(MAE(test_observations, model1_test_pred),
                                   MAE(test_observations, model2_test_pred)),
                     "trainR2" = c(R2(train_observations, model1_train_pred),
                                   R2(train_observations, model2_train_pred)),
                     "testR2" = c(R2(test_observations, model1_test_pred),
                                  R2(test_observations, model2_test_pred))
  )
  )
  
}


# Preparing Model Variables -----------------------------------------------

metadata <- c("UID", "Gene", "Guide", "LFC")

predictors <- setdiff(names(variables), c(metadata, "centered_LFC"))

# Training Linear Model ---------------------------------------------------

# Preparing formula for modelling
lm_formula <- as.formula(paste( "centered_LFC ~ ",
                                paste(predictors, collapse = " + ")))

linear_model <- lm(lm_formula ,
                   data = variables[train_index, ])


# Evaluating Model --------------------------------------------------------

summary(linear_model)

# Residual standard error: 0.3255 on 41459 degrees of freedom
# Multiple R-squared:  0.04302,	Adjusted R-squared:  0.04152 
# F-statistic: 28.67 on 65 and 41459 DF,  p-value: < 2.2e-16


# Adding Interaction Terms Between Positions ------------------------------

adjacent_terms <- paste0(
  "(pos", 1:19, "_A + pos", 1:19, "_C + pos", 1:19, "_G + pos", 1:19, "_T) * ",
  "(pos", 2:20, "_A + pos", 2:20, "_C + pos", 2:20, "_G + pos", 2:20, "_T)"
)

# Preparing formula for modelling
interaction_lm_formula <- as.formula(paste(
  "centered_LFC ~",
  paste(predictors, collapse = " + "),
  "+",
  paste(adjacent_terms, collapse = " + ")
)
)

# Training new model with cross-terms
lm_interaction_terms <- lm(interaction_lm_formula ,
                           data = variables[train_index, ]
)

summary(lm_interaction_terms)

# Residual standard error: 0.3236 on 41288 degrees of freedom
# Multiple R-squared:  0.05838,	Adjusted R-squared:  0.053 
# F-statistic: 10.85 on 236 and 41288 DF,  p-value: < 2.2e-16


# Comparing models --------------------------------------------------------

comparison <- get_summary(linear_model, lm_interaction_terms, 
                          "linear_model", "linear_model_w_cross_terms")

# Saving Linear Model -----------------------------------------------------

saveRDS(lm_interaction_terms, 
        "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")
