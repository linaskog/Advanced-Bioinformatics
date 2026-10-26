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

get_summary <- function(models, names) {
  
  train_observations <- variables[train_index, "centered_LFC"]
  test_observations <- variables[test_index, "centered_LFC"]
  
  df <- data.frame(Model = names)
  
  for (i in 1:length(names)) {
    
    train_pred <- predict(models[[i]], variables[train_index,])
    test_pred <- predict(models[[i]], variables[test_index, ])
    
    df$R2_train[i] <- basicstats::R2(train_observations, train_pred)
    df$R2_test[i] <- basicstats::R2(test_observations, test_pred)
    df$RMSE_train[i] <- basicstats::RMSE(train_observations, train_pred)
    df$RMSE_test[i] <- basicstats::RMSE(test_observations, test_pred)
    df$MAE_train[i] <- basicstats::MAE(train_observations, train_pred)
    df$MAE_test[i] <- basicstats::MAE(test_observations, test_pred)
    df$Spearman_train[i] <- basicstats::correlation(train_observations, 
                                                    train_pred,
                                                    "spearman")
    df$Spearman_test[i] <- basicstats::correlation(test_observations,
                                                   test_pred,
                                                   "spearman")
  }
  
  return( df )
  }


# Preparing Model Variables -----------------------------------------------

metadata <- c("UID", "Gene", "Guide", "LFC")

predictors <- setdiff(names(variables), c(metadata, "centered_LFC"))


# Creating Baseline Model for Comparison ----------------------------------

baseline <- lm(centered_LFC ~ GC,
               data = variables[train_index, ])

summary(baseline)

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

# Adding Additional Interaction Terms -------------------------------------

triple_adj_terms <- paste0(
  "(pos", 1:18, "_A + pos", 1:18, "_C + pos", 1:18, "_G + pos", 1:18, "_T) * ",
  "(pos", 2:19, "_A + pos", 2:19, "_C + pos", 2:19, "_G + pos", 2:19, "_T) * ",
  "(pos", 3:20, "_A + pos", 3:20, "_C + pos", 3:20, "_G + pos", 3:20, "_T)"
  )

# Preparing formula for modelling additional interactions
add_inter_lm_formula <- as.formula(paste(
  "centered_LFC ~ ",
  paste(predictors, collapse = " + "),
  "+",
  paste(adjacent_terms, collapse = " + "),
  "+",
  paste(triple_adj_terms, collapse = " + ")
  ))

# Training new model with additional cross-terms
add_inter_lm <- lm(add_inter_lm_formula,
                   data = variables[train_index,])

summary(add_inter_lm)

# Residual standard error: 0.3222 on 40640 degrees of freedom
# Multiple R-squared:  0.08071,	Adjusted R-squared:  0.06071 
# F-statistic: 4.036 on 884 and 40640 DF,  p-value: < 2.2e-16

# Comparing models --------------------------------------------------------

comparison <- get_summary(list(baseline,
                               linear_model, 
                            lm_interaction_terms,
                            add_inter_lm),
                          c("Baseline (only GC)",
                            "Linear Model (No cross terms",
                            "LM Dinucleotide Cross-Terms",
                            "LM Trinucleotide Cross-Terms")
                          )

comparison

# Model                         R2_train    R2_test     RMSE_train RMSE_test 
# Baseline (only GC)            0.001837295 0.001798662  0.3321744 0.3344780 
# Linear Model (No cross terms  0.04302252  0.04121038  0.3252493  0.3278084 
# LM Dinucleotide Cross-Terms   0.05837949  0.04644043  0.3226290  0.3269131 
# LM Trinucleotide Cross-Terms  0.08070696  0.03678290  0.3187810  0.3285644 

# Model                         MAE_train  MAE_test Spearman_train Spearman_test
# Baseline (only GC)            0.2411301 0.2417564 0.04075234     0.04586712
# Linear Model (No cross terms  0.2375401 0.2383944 0.2016518      0.2012359
# LM Dinucleotide Cross-Terms   0.2360659 0.2380683 0.2325099      0.2163131
# LM Trinucleotide Cross-Terms  0.2342924 0.2402606 0.2666679      0.2043376


# Saving Linear Model -----------------------------------------------------

saveRDS(lm_interaction_terms, 
        "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")
