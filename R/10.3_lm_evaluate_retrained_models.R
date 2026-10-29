# 10.3_lm_evaluate_retrained_models.R

# Purpose:
# Perform additional evaluation of linear models produced in R/10_linear_model.R

# Input:
# Sources "R10.2_lm_sign_coeff_train.R"


# Loading Data and Libraries ----------------------------------------------

source("R/10.2_lm_sign_coeff_train.R")


# Function for Comparing Models -------------------------------------------

get_summary <- function(models, names) {
  
  # Extract observed centered LFC data
  train_observations <- variables[train_index, "centered_LFC"]
  test_observations <- variables[test_index, "centered_LFC"]
  
  df <- data.frame(Model = names)
  
  for (i in 1:length(names)) {
    
    # Predict new centrered LFC data
    train_pred <- predict(models[[i]], variables[train_index,])
    test_pred <- predict(models[[i]], variables[test_index, ])
    
    # Compute model evaluation metrics
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


# Comparing Models --------------------------------------------------------

retrained_lms <- get_summary(list(sig_LM,
                                  sig_LM_di,
                                  sig_LM_tri),
                             c("Linear Model",
                               "LM + Pair Cross-Terms",
                               "LM + Pair + Triplet Cross-Terms"))

retrained_lms
#                             Model   R2_train    R2_test RMSE_train RMSE_test
#                    Linear Model 0.04243236 0.04081382  0.3253495 0.3278762
#           LM + Pair Cross-Terms 0.05326290 0.04509475  0.3235044 0.3271437
# LM + Pair + Triplet Cross-Terms 0.04952743 0.03827260  0.3241420 0.3283102

#                         Model MAE_train  MAE_test Spearman_train Spearman_test
#                    Linear Model 0.2375673 0.2384014   0.2005482     0.2002227
#           LM + Pair Cross-Terms 0.2365310 0.2382005   0.2220417     0.2100304
# LM + Pair + Triplet Cross-Terms 0.2369924 0.2384865   0.2104659     0.1959135

# Exploring New  Coefficients ----------------------------------------------

# Extracting regression coefficients
LM_coeffs <- summary(sig_LM)$coefficients
LM_di_coeffs <- summary(sig_LM_di)$coefficients
LM_tri_coeffs <- summary(sig_LM_tri)$coefficients


# Extracting top 10 regression coefficients by Pr(>|t|)
top_10_coeffs <- data.frame(
  LM = rownames(LM_coeffs[order(LM_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_di = rownames(LM_di_coeffs[order(LM_di_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_tri = rownames(LM_tri_coeffs[order(LM_tri_coeffs[,"Pr(>|t|)"]), ])[1:10]
  
)

# Extracting bottom 10 regression coefficients by Pr(>|t|)
bottom_10_coeffs <- data.frame(
  LM = rownames(LM_coeffs[order(-LM_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_di = rownames(LM_di_coeffs[order(-LM_di_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_tri = rownames(LM_tri_coeffs[order(-LM_tri_coeffs[,"Pr(>|t|)"]), ])[1:10]
  
)
