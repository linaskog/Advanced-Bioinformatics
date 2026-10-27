# 13.1_xgb_model_variables.R

# Purpose:
# Investigate xgb model variables

# Input:
# - sources data from "R/13_xgb_model_evaluation.R"


# Reading in Libraries ----------------------------------------------------

library(xgboost)
library(SHAPforxgboost)
#install.packages("R/basicstats", repos = NULL, type = "source")
library(basicstats)

# Reading in Data ---------------------------------------------------------

xgb_model <- xgb.load("~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

y <- as.numeric(variables$centered_LFC)

pos <- as.matrix(variables[, grep("pos", colnames(variables))])
additional_data <- as.matrix(variables[, c("GC", "strand", "chromosome", "region_type",
                                           "gene_overlap_count", "H3K27ac_peak_overlap")
])

X <- cbind(pos, additional_data)

train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds"
)

# SHAP analysis ----------------------------------------------------------

shap_values <- shap.values(
  xgb_model = xgb_model,
  X_train = X[train_index, ]
)

shap_long <- shap.prep(shap_contrib = shap_values$shap_score, 
                       X_train = X[train_index,])

shap.plot.summary(shap_long, dilute = 100)


# Importance --------------------------------------------------------------


# Saving Final XGB Model --------------------------------------------------

xgb.save(xgb_model, "~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")

