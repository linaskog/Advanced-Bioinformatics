# 13.2_xgb_model_genes.R

# Purpose:
# Investigate genes favored by sgRNA score in xgb_model


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


# Computing Table For training sgRNA Efficacy --------------------------------

positions <- as.matrix(variables[, grep("pos", colnames(variables))])
additional_data <- as.matrix(
  variables[, c("GC", "strand", "chromosome", "region_type",
                "gene_overlap_count", "H3K27ac_peak_overlap")
            ]
  )

X <- cbind(positions, additional_data)
rownames(X) <- variables$UID


training <- variables[train_index, 
                      c("UID", "Gene", "centered_LFC")]
training$predictions <- predict(xgb_model, 
                                X[train_index, ])

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


# Printing Genes Of Interest ----------------------------------------------

good_train_obs <- gene_avg_score[order(gene_avg_score$Observed), "Gene"][1:10]
print(paste("Guides in training set with high observed actvity:", 
            paste(good_train_obs, collapse = ", ")))
poor_train_obs <- gene_avg_score[order(-gene_avg_score$Observed), "Gene"][1:10]
print(paste("Guides in training set with low observed actvity:", 
            paste(poor_train_obs, collapse = ", ")))

good_train_pred <- gene_avg_score[order(gene_avg_score$Predicted), "Gene"][1:10]
print(paste("Guides in training set with high predicted actvity:", 
            paste(good_train_pred, collapse = ", ")))
poor_train_pred <- gene_avg_score[order(-gene_avg_score$Predicted), "Gene"][1:10]
print(paste("Guides in training set with low predicted actvity:", 
            paste(poor_train_pred, collapse = ", ")))

low_res_train <- gene_avg_score[order(gene_avg_score$Residual), "Gene"][1:10]
print(paste("Guides in training set that were accurately predicted:", 
            paste(low_res_train, collapse = ", ")))
high_res_train <- gene_avg_score[order(-gene_avg_score$Residual), "Gene"][1:10]
print(paste("Guides in training set that were not accurately predicted:", 
            paste(high_res_train, collapse = ", ")))

# Computing Table For Testing sgRNA Efficacy --------------------------------

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


# Printing Genes Of Interest ----------------------------------------------

good_test_obs <- gene_avg_score_test[order(gene_avg_score_test$Observed), 
                                     "Gene"][1:10]
print(paste("Guides in testing set with high observed actvity:", 
            paste(good_test_obs, collapse = ", ")))
poor_test_obs <- gene_avg_score_test[order(-gene_avg_score_test$Observed), 
                                     "Gene"][1:10]
print(paste("Guides in testing set with low observed actvity:", 
            paste(poor_test_obs, collapse = ", ")))

good_test_pred <- gene_avg_score_test[order(gene_avg_score_test$Predicted), 
                                      "Gene"][1:10]
print(paste("Guides in testing set with high predicted actvity:", 
            paste(good_test_pred, collapse = ", ")))
poor_test_pred <- gene_avg_score_test[order(-gene_avg_score_test$Predicted), 
                                      "Gene"][1:10]
print(paste("Guides in testing set with low predicted actvity:", 
            paste(poor_test_pred, collapse = ", ")))

low_res_test <- gene_avg_score_test[order(gene_avg_score_test$Residual), 
                                    "Gene"][1:10]
print(paste("Guides in testing set that were accurately predicted:", 
            paste(low_res_test, collapse = ", ")))
high_res_test <- gene_avg_score_test[order(-gene_avg_score_test$Residual), 
                                     "Gene"][1:10]
print(paste("Guides in testing set that were not accurately predicted:", 
            paste(high_res_test, collapse = ", ")))


