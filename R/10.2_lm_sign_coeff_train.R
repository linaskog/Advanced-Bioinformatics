# 10.2_lm_sign_coeff_train.R

# Purpose:
# Explore Coefficients For the Linear Models from R/10_linear_model.R

# Input:
# Sources "R/10_linear_model.R"

# Output:

# Loading Data and Libraries ----------------------------------------------

source("R/10_linear_model.R")

# Exploring Coefficients --------------------------------------------------

LM_coeffs <- summary(linear_model)$coefficients
LM_di_coeffs <- summary(lm_interaction_terms)$coefficients
LM_tri_coeffs <- summary(add_inter_lm)$coefficients

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

# Finding Significant Coefficients ----------------------------------------

cut_off <- 0.1

# Finding significant coefficinets for linear model
nrow(LM_coeffs)
# 66 variables originally
sig_LM_coeffs <- rownames(LM_coeffs[LM_coeffs[, 'Pr(>|t|)'] < sig_level, ])
sig_LM_coeffs <- setdiff(sig_LM_coeffs, "(Intercept)")
length(sig_LM_coeffs) 
# Down to 35 significant variables

# Finding significant coefficinets for linear model + pair cross-terms
nrow(LM_di_coeffs)
# 237 variables originally
sig_LM_di_coeffs <- rownames(LM_di_coeffs[
  LM_di_coeffs[ ,'Pr(>|t|)'] < sig_level, ])
sig_LM_di_coeffs <- setdiff(sig_LM_di_coeffs, "(Intercept)")
length(sig_LM_di_coeffs)
# Down to 94 significant variables

# Finding significant coefficinets for linear model + pair + triplet cross-terms
nrow(LM_tri_coeffs)
# 885 variables originally
sig_LM_tri_coeffs <- rownames(LM_tri_coeffs[
  LM_tri_coeffs[ ,'Pr(>|t|)'] < sig_level, ])
sig_LM_tri_coeffs <- setdiff(sig_LM_tri_coeffs, "(Intercept)")
length(sig_LM_tri_coeffs)
# Down to 135 significant variables

# Training Linear Models Using Only Significant Variables -----------------

# Linear Model 
sig_LM_formula <- as.formula(
  paste0(
    "centered_LFC ~ ",
    paste0(
      sig_LM_coeffs, 
      collapse = " + "
    )
  )
)

sig_LM <- lm(sig_LM_formula,
             data = variables[train_index,])

summary(sig_LM)

# Residual standard error: 0.3255 on 41489 degrees of freedom
# Multiple R-squared:  0.04243,	Adjusted R-squared:  0.04162 
# F-statistic: 52.53 on 35 and 41489 DF,  p-value: < 2.2e-16


# Linear Model + Dinucleotide Cross-Terms
sig_LM_di_formula <- as.formula(
  paste0(
    "centered_LFC ~ ",
    paste0(
      sig_LM_di_coeffs, 
      collapse = " + "
    )
  )
)

sig_LM_di <- lm(sig_LM_di_formula,
                data = variables[train_index,])

summary(sig_LM_di)

# Residual standard error: 0.3239 on 41430 degrees of freedom
# Multiple R-squared:  0.05326,	Adjusted R-squared:  0.05111 
# F-statistic:  24.8 on 94 and 41430 DF,  p-value: < 2.2e-16

# Linear Model + Dinucleotide + Trinucleotide Cross-Terms
sig_LM_tri_formula <- as.formula(
  paste0(
    "centered_LFC ~ ",
    paste0(
      sig_LM_tri_coeffs, 
      collapse = " + "
    )
  )
)

sig_LM_tri <- lm(sig_LM_tri_formula,
                 data = variables[train_index,])

summary(sig_LM_tri)

# Residual standard error: 0.3247 on 41389 degrees of freedom
# Multiple R-squared:  0.04953,	Adjusted R-squared:  0.04643 
# F-statistic: 15.98 on 135 and 41389 DF,  p-value: < 2.2e-16