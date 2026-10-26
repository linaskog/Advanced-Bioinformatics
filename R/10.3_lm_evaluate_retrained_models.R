# 10.3_lm_evaluate_retrained_models.R

# Purpose:
# Perform additional evaluation of linear models produced in R/10_linear_model.R

# Input:
# Sources "R10.2_lm_sign_coeff_train.R"


# Loading Data and Libraries ----------------------------------------------

source("R/10.2_lm_sign_coeff_train.R")


# Function for Comparing Models -------------------------------------------

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


# Comparing Models --------------------------------------------------------

retrained_lms <- get_summary(list(sig_LM,
                                  sig_LM_di,
                                  sig_LM_tri),
                             c("Linear Model",
                               "LM + Dinucleotide Cross-Terms",
                               "LM + Di + Trinucleotide Cross-Terms"))

# Exploring New  Coefficients ----------------------------------------------

LM_coeffs <- summary(sig_LM)$coefficients
LM_di_coeffs <- summary(sig_LM_di)$coefficients
LM_tri_coeffs <- summary(sig_LM_tri)$coefficients

top_10_coeffs <- data.frame(
  LM = rownames(LM_coeffs[order(LM_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_di = rownames(LM_di_coeffs[order(LM_di_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_tri = rownames(LM_tri_coeffs[order(LM_tri_coeffs[,"Pr(>|t|)"]), ])[1:10]
  
)

bottom_10_coeffs <- data.frame(
  LM = rownames(LM_coeffs[order(-LM_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_di = rownames(LM_di_coeffs[order(-LM_di_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_tri = rownames(LM_tri_coeffs[order(-LM_tri_coeffs[,"Pr(>|t|)"]), ])[1:10]
  
)
