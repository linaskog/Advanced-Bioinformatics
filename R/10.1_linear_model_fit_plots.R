# 10.1_linear_model_fit_plots.R

# Purpose:
# Perform additional evaluation of linear models produced in R/10_linear_model.R

# Input:
# Sources "R/10_linear_model.R"


# Loading Data and Libraries ----------------------------------------------

source("R/10_linear_model.R")
library(ggplot2)

# Preparing data frame for observations and predictions of training data set
train_data <- data.frame(Observations = variables[train_index, "centered_LFC"],
                         LM_predictions = predict(linear_model,
                                                  variables[train_index,]),
                         di_LM_predictions = predict(lm_interaction_terms,
                                                     variables[train_index,]),
                         tri_LM_predictions = predict(add_inter_lm,
                                                      variables[train_index,]))

# Preparing data frame for observations and predictions of test data set
test_data <- data.frame(Observations = variables[test_index, "centered_LFC"],
                         LM_predictions = predict(linear_model,
                                                  variables[test_index,]),
                         di_LM_predictions = predict(lm_interaction_terms,
                                                     variables[test_index,]),
                         tri_LM_predictions = predict(add_inter_lm,
                                                      variables[test_index,]))

# Plotting Fit of Linear Models on  Training Data -------------------------

# Normal Linear Model fit
ggplot(data = train_data, 
       aes(x = LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM Train Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "blue") + 
  theme_bw()

# Linear model + pair cross terms
ggplot(data = train_data, 
       aes(x = di_LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM + Pair Cross-Terms Train Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "blue") + 
  theme_bw()

# Linear model + pair + triplet cross terms
ggplot(data = train_data, 
       aes(x = tri_LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM + Triplet Cross-Terms Train Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "blue") + 
  theme_bw()


# Plotting Fit of Linear Models on  Testing Data -------------------------

# Normal Linear Model fit
ggplot(data = test_data, 
       aes(x = LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM Test Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw()

# Linear model + pair cross terms
ggplot(data = test_data, 
       aes(x = di_LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM + Pair Cross-Terms Test Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw()

# Linear model + pair + triplet cross terms
ggplot(data = test_data, 
       aes(x = tri_LM_predictions,
           y = Observations)) +
  geom_point() +
  labs(title = "LM + Triplet Cross-Terms Test Pred vs Obs") + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw()
