# 10.2_lm_sign_coeff_train.R

# Purpose:
# Attempting to retrain linear models only on significant coefficients

# Input:
# Sources "R/10_linear_model.R"

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

top_10_coeffs
# LM        LM + pairs        LM + pairs + triplets
#	pos20_G	  pos19_G	          pos20_G
#	pos18_C	  pos20_G	          pos19_G
#	pos19_G	  pos19_G:pos20_C	  pos19_G:pos20_C
#	pos17_G	  pos19_C	          (Intercept)
#	pos19_C	  pos18_C	          H3K27ac_peak_overlap
#	pos17_C	  (Intercept)	      pos18_C:pos20_G
#	pos17_A	  pos17_G	          pos18_C
#	pos19_A	  pos17_C	          pos19_C:pos20_G
#	pos20_C	  pos18_G	          pos1_C:pos2_C
#	pos10_A	  pos18_A	          pos19_C

bottom_10_coeffs <- data.frame(
  LM = rownames(LM_coeffs[order(-LM_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_di = rownames(LM_di_coeffs[order(-LM_di_coeffs[,"Pr(>|t|)"]), ])[1:10],
  LM_tri = rownames(LM_tri_coeffs[order(-LM_tri_coeffs[,"Pr(>|t|)"]), ])[1:10]
  
)
bottom_10_coeffs
#                   LM              LM_di                  LM_tri
#              pos1_A      pos6_C:pos7_A          pos9_A:pos10_A
#              pos7_C      pos1_C:pos2_C           pos7_A:pos8_G
#  gene_overlap_count      pos5_G:pos6_C         pos11_G:pos13_C
#              pos8_C gene_overlap_count pos18_C:pos19_G:pos20_A
#             pos16_G    pos10_G:pos11_A           pos4_G:pos6_C
#              pos3_A      pos5_G:pos6_G           pos5_C:pos7_C
#              pos2_C      pos3_A:pos4_C    pos1_G:pos2_G:pos3_A
#              pos9_C    pos16_C:pos17_G           pos3_A:pos4_C
#              pos6_A    pos16_A:pos17_A pos11_G:pos12_C:pos13_A
#        region_type      pos4_C:pos5_A   pos8_C:pos9_G:pos10_G

# Finding Significant Coefficients ----------------------------------------

sig_level <- 0.1

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