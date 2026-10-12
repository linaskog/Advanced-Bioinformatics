# 14_final_model_exploration.R

# Purpose:
# Explore final linear and xgb models

# Input:
# - Linear model
# - XGB model
# - variables for training, training index and testing index


# Reading in necessary libraries ------------------------------------------

library(xgboost)
library(ggplot2)

# Reading in Models -------------------------------------------------------

xgb_model <- xgb.load("~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")
linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")


# Reading in data ---------------------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")
train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")


# Visualizing Model Fit ---------------------------------------------------


