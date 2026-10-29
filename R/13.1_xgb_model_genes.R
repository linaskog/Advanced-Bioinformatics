# 13.1_xgb_model_genes.R

# Purpose:
# Investigating observed and model predicted differences for genes marked as
# essential

# Input:
# - variables
# - training and testing indices used for training model
# - trained xgb model
# - essential genes computed in 10.4_linear_model_genes.R

# Output:
# - distribution plots of gene average predicted and observed centered LFC 
#   scores and residuals "results/models/"

# Import variables and models ---------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")

xgb_model <- xgb.load(
  "~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj"
)


# Prepare Modelling Data --------------------------------------------------


essential_genes <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/essential_genes.rds"
)

pos <- as.matrix(variables[, grep("pos", colnames(variables))])
additional_data <- as.matrix(variables[, c("GC", "strand", "chromosome", "region_type",
                                           "gene_overlap_count", "H3K27ac_peak_overlap")
])

X <- cbind(pos, additional_data)


# Investigating Training Dataset ------------------------------------------

# Extracting guides from training data set
training <- variables[train_index, 
                      c("UID", "Gene", "centered_LFC")]

training$predictions <- predict(xgb_model, 
                                X[train_index, ])

# Computing absolute residual vale between predicted and observed centered LFC
training$residual <- abs(training$predictions - training$centered_LFC)

# Computing Gene Averaged Metrics -----------------------------------------

# Average observed score per gene
train_obs <- tapply(training$centered_LFC,
                    training$Gene, 
                    mean,
                    na.rm = TRUE)

# Average predicted score per gene
train_pred <- tapply(training$predictions,
                     training$Gene, 
                     mean,
                     na.rm = TRUE)

# Average residual score per gene
train_res <- tapply(training$residual,
                    training$Gene, 
                    mean,
                    na.rm = TRUE)

# Final Data Frame With All Metrics
train_genes <- unique(training$Gene)
gene_avg_score_train <- data.frame(
  "Gene" = train_genes,
  "Observed" = train_obs[train_genes],
  "Predicted" = train_pred[train_genes],
  "Residual" = train_res[train_genes]
)

# Adding column for gene essentiality
gene_avg_score_train$Essential <- ifelse(gene_avg_score_train$Gene %in% essential_genes$Gene,
                                   TRUE,
                                   FALSE)

# Visualizing Training Set Gene Score Distributions -----------------------

# Plotting average residual score distribution for train set genes
essential_gene_res_tr <- ggplot(gene_avg_score_train, 
                                aes(x = Residual, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Train Residual Score Distribution",
       x = "Residual Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_res_tr.png",
        plot = essential_gene_res_tr,
        width = 10, height = 6, dpi = 300)

# Plotting average predicted score distribution for train set genes
essential_gene_pred_tr <- ggplot(gene_avg_score_train, 
                                 aes(x = Predicted, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Train Predicted Centered LFC Score Distribution",
       x = "Predicted Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_pred_tr.png",
        plot = essential_gene_pred_tr,
        width = 10, height = 6, dpi = 300)

# Plotting average observed score distribution for train set genes
essential_gene_obs_tr <- ggplot(gene_avg_score_train, 
                                aes(x = Observed, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Train Observed Centered LFC Score Distribution",
       x = "Observed Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_obs_tr.png",
        plot = essential_gene_obs_tr,
        width = 10, height = 6, dpi = 300)



# Investigating Training Dataset ------------------------------------------

# Extracting guides from training data set
testing <- variables[test_index, 
                     c("UID", "Gene", "centered_LFC")
]
testing$predictions <- predict(xgb_model,
                               X[test_index, ])

testing$residual <- abs(testing$predictions - testing$centered_LFC)


# Computing Gene Averaged Metrics -----------------------------------------

# Average observed score per gene
test_obs <- tapply(testing$centered_LFC,
                   testing$Gene,
                   mean,
                   na.rm = TRUE)

# Average predicted score per gene
test_pred <- tapply(testing$predictions,
                    testing$Gene,
                    mean,
                    na.rm = TRUE)

# Average residual score per gene
test_res <- tapply(testing$residual,
                   testing$Gene,
                   mean,
                   na.rm = TRUE)

# Final Data Frame With All Metrics
test_genes <- unique(testing$Gene)
gene_avg_score_test <- data.frame(
  "Gene" = test_genes,
  "Observed" = test_obs[test_genes],
  "Predicted" = test_pred[test_genes],
  "Residual" = test_res[test_genes]
)

gene_avg_score_test$Essential <- ifelse(
  gene_avg_score_test$Gene %in% essential_genes$Gene,
  TRUE,
  FALSE)


# Visualizing Training Set Gene Score Distributions -----------------------

# Plotting average residual score distribution for test set genes
essential_gene_res_te <- ggplot(gene_avg_score_test, 
                                aes(x = Residual, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Test Residual Centered LFC Score Distribution",
       x = "Residual Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_res_te.png",
        plot = essential_gene_res_te,
        width = 10, height = 6, dpi = 300)

# Plotting average predicted score distribution for test set genes
essential_gene_pred_te <-ggplot(gene_avg_score_test, 
                                aes(x = Predicted, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Test Predicted Centered LFC Score Distribution",
       x = "Predicted Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_pred_te.png",
        plot = essential_gene_pred_te,
        width = 10, height = 6, dpi = 300)

# Plotting average observed score distribution for test set genes
essential_gene_obs_te <- ggplot(gene_avg_score_test, 
                                aes(x = Observed, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "XGB Test Observed Centered LFC Score Distribution",
       x = "Observed Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/xgb_ess_gene_obs_te.png",
        plot = essential_gene_obs_te,
        width = 10, height = 6, dpi = 300)



# Function for computing metrics ------------------------------------------

compute_metrics <- function(data, name) {
  
  df <- data.frame(
    "Dataset" = name,
    "Minimum" = min(data),
    "Maximum" = max(data),
    "Mean" = mean(data),
    "StandardDeviation" = sd(data)
  )
}

# Residual Metrics
train_res_essential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == TRUE, "Residual"], "TrainingEssential")
train_res_nonessential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == FALSE, "Residual"], "TrainingNonEssential")
test_res_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Residual"], "TestEssential")
test_res_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Residual"], "TestgNonEssential")

res_metrics <- rbind(train_res_essential,
                     train_res_nonessential,
                     test_res_essential,
                     test_res_nonessential)
res_metrics
#                Dataset      Minimum  Maximum      Mean StandardDeviation
#    TrainingEssential 0.0016442754 2.199489 0.3419923         0.2475622
# TrainingNonEssential 0.0008247302 1.958313 0.1951682         0.1200298
#        TestEssential 0.0008239819 2.068420 0.3451042         0.2642991
#    TestgNonEssential 0.0047005108 1.370204 0.1974192         0.1175499

# Predicted Metrics
train_pred_essential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == TRUE, "Predicted"], "TrainingEssential")
train_pred_nonessential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == FALSE, "Predicted"], "TrainingNonEssential")
test_pred_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Predicted"], "TestEssential")
test_pred_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Predicted"], "TestNonEssential")

pred_metrics <- rbind(train_pred_essential,
                      train_pred_nonessential,
                      test_pred_essential,
                      test_pred_nonessential)
pred_metrics
#               Dataset    Minimum   Maximum         Mean StandardDeviation
#    TrainingEssential -0.2318455 0.1735972 -0.014058247        0.04134573
# TrainingNonEssential -0.1731526 0.2238040  0.004562518        0.04084829
#        TestEssential -0.1484434 0.1133934 -0.016010231        0.04004724
#     TestNonEssential -0.2014117 0.2121579  0.005881400        0.04282813

# Observed Metrics
train_obs_essential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == TRUE, "Observed"], "TrainingEssential")
train_obs_nonessential <- compute_metrics(gene_avg_score_train[
  gene_avg_score_train$Essential == FALSE, "Observed"], "TrainingNonEssential")
test_obs_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Observed"], "TestEssential")
test_obs_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Observed"], "TestNonEssential")

obs_metrics <- rbind(train_obs_essential,
                     train_obs_nonessential,
                     test_obs_essential,
                     test_obs_nonessential)
obs_metrics
#                Dataset       Minimum      Maximum          Mean StandardDeviation
#    TrainingEssential -1.480659e-16 2.220446e-16 -1.834013e-19      2.880866e-17
# TrainingNonEssential -1.110223e-16 7.403294e-17  8.354851e-20      1.196423e-17
#        TestEssential -1.480297e-16 2.220446e-16 -1.020386e-18      3.215154e-17
#     TestNonEssential -7.401939e-17 5.551115e-17  9.751321e-20      1.153447e-17
