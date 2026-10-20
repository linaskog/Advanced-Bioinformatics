# 11_linear_model_evaluation.R

# Purpose:
# Perform additional evaluation of linear model produced in 
# R/10_linear_model.R

# Input:
# - linear model as RDS object
# - variables as RDS object
# - train and test indices as RDS objects

# Reading in Libraries ----------------------------------------------------

library(ggplot2)
library(caret)

# Reading in data and models ----------------------------------------------

linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")

test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")


# Producing Diagnostic Plots ----------------------------------------------

par(mfrow = c(2, 2))
plot(linear_model)

# Visualizing Predictive Performance --------------------------------------

train_data <- data.frame(obs = variables[train_index, "centered_LFC"])
train_data$pred <- predict(linear_model, variables[train_index,])

test_data <- data.frame(obs = variables[test_index, "centered_LFC"])
test_data$pred <- predict(linear_model, variables[test_index,])

ggplot(data = train_data, aes(obs, pred)) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(title = "Training Predicted vs. Observed",
       x = "Observed",
       y = "Predicted") +
  theme_bw()

ggplot(data = test_data, aes(obs, pred)) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(title = "Test Predicted vs. Observed",
       x = "Observed",
       y = "Predicted") +
  theme_bw()



# Finding Important Variables ---------------------------------------------

coeff_summary <- summary(linear_model)$coefficients

coeff_summary <- coeff_summary[order(coeff_summary[, 3]), ]

coeff_summary <- coeff_summary[setdiff(rownames(coeff_summary), 
                                             "(Intercept)"), ]
top_ten <- coeff_summary[1:10,]
