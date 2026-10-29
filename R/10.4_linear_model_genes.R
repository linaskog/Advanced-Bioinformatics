# 10.4_linear_model_genes.R

# Purpose:
# Investigating observed and model predicted differences for genes marked as
# essential

# Input:
# - variables
# - training and testing indices used for training model
# - trained linear model

# Output:
# - RDS object of essential genes with columns gene name and gene average LFC
# - distribution plots of gene average predicted and observed centered LFC 
#   scores and residuals "results/models/"

# Import variables and models ---------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")

train_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/train_index.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")

linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds"
)


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

# Computing Essential Genes -----------------------------------------------

# Computing gene averaged Log2 Fold Change
gene_averaged_LFC <- tapply(variables$LFC,
                            variables$Gene,
                            mean,
                            na.rm = TRUE)

# Initializing data frame for storing essential genes
essential_genes <- data.frame("Gene" = unique(variables$Gene)) 
essential_genes$avg_LFC <- gene_averaged_LFC[essential_genes$Gene]
# total genes: 17969
nrow(essential_genes)

# Selecting "essential" genes as those scoring in the bottom 25% of LFC
q25 <- quantile(gene_averaged_LFC, 0.25)
essential_genes <- essential_genes[essential_genes$avg_LFC <= q25, ]

# Computing essential gene metrics
# Number essential genes: 4493
nrow(essential_genes)

min(essential_genes$avg_LFC)
mean(essential_genes$avg_LFC)
max(essential_genes$avg_LFC)
# min:-3.03954, mean:-0.459607, max:-0.1412059

# saving essential gene data frame for xgb model exploration
saveRDS(essential_genes, "data/processed/essential_genes.rds")


# Investigating Training Dataset ------------------------------------------

# Extracting guides from training data set
training <- variables[train_index, 
                      c("UID", "Gene", "centered_LFC")]

training$predictions <- predict(linear_model, 
                                variables[train_index, ])

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
gene_avg_score <- data.frame(
  "Gene" = train_genes,
  "Observed" = train_obs[train_genes],
  "Predicted" = train_pred[train_genes],
  "Residual" = train_res[train_genes]
)

# Adding column for gene essentiality
gene_avg_score$Essential <- ifelse(gene_avg_score$Gene %in% essential_genes$Gene,
                                   TRUE,
                                   FALSE)


# Analysing gene averaged metrics -----------------------------------------

# Total number genes in training data set: 14375 genes
nrow(gene_avg_score)

# Number essential genes in training data set: 3600
nrow(gene_avg_score[gene_avg_score$Essential == TRUE,])
# Fraction essential genes in training data set: 0.25
round(nrow(gene_avg_score[gene_avg_score$Essential == TRUE,]
           ) / nrow(gene_avg_score), 2)

# Number essential genes in training data set: 10775
nrow(gene_avg_score[gene_avg_score$Essential == FALSE,])
# Fraction essential genes in training data set: 0.75
round(nrow(gene_avg_score[gene_avg_score$Essential == FALSE,] 
           ) / nrow(gene_avg_score), 2)


# Visualizing Training Set Gene Score Distributions -----------------------

# Plotting average residual score distribution for train set genes
essential_gene_res_tr <- ggplot(gene_avg_score, 
                                aes(x = Residual, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Training Residual Score Distribution",
       x = "Residual Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_res_tr.png",
        plot = essential_gene_res_tr,
        width = 10, height = 6, dpi = 300)

# Plotting average predicted score distribution for train set genes
essential_gene_pred_tr <- ggplot(gene_avg_score, 
                                 aes(x = Predicted, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Training Predicted Centered LFC Score Distribution",
       x = "Predicted Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_pred_tr.png",
        plot = essential_gene_pred_tr,
        width = 10, height = 6, dpi = 300)

# Plotting average observed score distribution for train set genes
essential_gene_obs_tr <- ggplot(gene_avg_score, aes(x = Observed, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Training Observed Centered LFC Score Distribution",
       x = "Observed Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_obs_tr.png",
        plot = essential_gene_obs_tr,
        width = 10, height = 6, dpi = 300)



# Investigating Training Dataset ------------------------------------------

# Extracting guides from training data set
testing <- variables[test_index, 
                     c("UID", "Gene", "centered_LFC")
                     ]
testing$predictions <- predict(linear_model,
                               variables[test_index, ])

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
# Total number genes in training data set: 3594
nrow(gene_avg_score_test)

# Number essential genes in training data set: 893
nrow(gene_avg_score_test[gene_avg_score_test$Essential == TRUE,])
# Fraction essential genes in training data set: 0.25
round(nrow(gene_avg_score_test[gene_avg_score_test$Essential == TRUE, ]
           ) / nrow(gene_avg_score_test), 2)

# Number essential genes in training data set: 2701
nrow(gene_avg_score_test[gene_avg_score_test$Essential == FALSE,])
# Fraction essential genes in training data set: 0.75
round(nrow(gene_avg_score_test[gene_avg_score_test$Essential == FALSE,]
           ) / nrow(gene_avg_score_test), 2)

# Visualizing Training Set Gene Score Distributions -----------------------

# Plotting average residual score distribution for test set genes
essential_gene_res_te <- ggplot(gene_avg_score_test, 
                                aes(x = Residual, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Test Residual Centered LFC Score Distribution",
       x = "Residual Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_res_te.png",
        plot = essential_gene_res_te,
        width = 10, height = 6, dpi = 300)

# Plotting average predicted score distribution for test set genes
essential_gene_pred_te <-ggplot(gene_avg_score_test, 
                               aes(x = Predicted, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Testing Predicted Centered LFC Score Distribution",
       x = "Predicted Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_pred_te.png",
        plot = essential_gene_pred_te,
        width = 10, height = 6, dpi = 300)

# Plotting average observed score distribution for test set genes
essential_gene_obs_te <- ggplot(gene_avg_score_test, 
                                aes(x = Observed, fill = Essential)) +
  geom_histogram() +
  theme_bw() + 
  labs(title = "LM Testing Observed Centered LFC Score Distribution",
       x = "Observed Centered LFC Score",
       y = "Count",
       fill = "Essential Gene")
ggsave( filename = "results/models/lm_ess_gene_obs_te.png",
        plot = essential_gene_obs_te,
        width = 10, height = 6, dpi = 300)



# Computing Metrics -------------------------------------------------------

# Residual Metrics
train_res_essential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == TRUE, "Residual"], "TrainingEssential")
train_res_nonessential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == FALSE, "Residual"], "TrainingNonEssential")
test_res_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Residual"], "TestEssential")
test_res_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Residual"], "TestgNonEssential")

res_metrics <- rbind(train_res_essential,
                     train_res_nonessential,
                     test_res_essential,
                     test_res_nonessential)

res_metrics

#               Dataset      Minimum  Maximum      Mean StandardDeviation
#    TrainingEssential 1.409265e-05 2.196987 0.3422662         0.2468379
# TrainingNonEssential 1.430359e-05 1.985004 0.1968571         0.1198260
#        TestEssential 7.631791e-04 2.043277 0.3448259         0.2623101
#    TestgNonEssential 1.826963e-03 1.387427 0.1988146         0.1171710

# Predicted Metrics
train_pred_essential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == TRUE, "Predicted"], "TrainingEssential")
train_pred_nonessential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == FALSE, "Predicted"], "TrainingNonEssential")
test_pred_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Predicted"], "TestEssential")
test_pred_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Predicted"], "TestNonEssential")

pred_metrics <- rbind(train_pred_essential,
                     train_pred_nonessential,
                     test_pred_essential,
                     test_pred_nonessential)
pred_metrics

# Observed Metrics
train_obs_essential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == TRUE, "Observed"], "TrainingEssential")
train_obs_nonessential <- compute_metrics(gene_avg_score[
  gene_avg_score$Essential == FALSE, "Observed"], "TrainingNonEssential")
test_obs_essential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == TRUE, "Observed"], "TestEssential")
test_obs_nonessential <- compute_metrics(gene_avg_score_test[
  gene_avg_score_test$Essential == FALSE, "Observed"], "TestNonEssential")

obs_metrics <- rbind(train_obs_essential,
                      train_obs_nonessential,
                      test_obs_essential,
                      test_obs_nonessential)
