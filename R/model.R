# 0?_model.R
# Advanced Bioinformatics: sgRNA quality project

# Purpose:
# Model sgRNA efficacy on sequence composition


# Reading in Packages -----------------------------------------------------
library(xgboost)
library(readxl)
library(tidyverse)
library(stringr)
library(DESeq2)


# Reading in Data ---------------------------------------------------------
excel_sheets("~/Advanced-Bioinformatics/41467_2018_5506_MOESM4_ESM.xlsx")

# Raw Count Data
gecko <- read_excel("~/Advanced-Bioinformatics/41467_2018_5506_MOESM4_ESM.xlsx", sheet = "C. GECKO_BC3_BJAB" )

# sgRNA Sequences
sequences <- read_excel("~/Advanced-Bioinformatics/41467_2018_5506_MOESM4_ESM.xlsx", sheet = "F. GeCKOlibrary2")

# Extracting Library B

lib_B_index <- grep("LibB", sequences$UID)

lib_B_sequences <- sequences[lib_B_index, ]

nrow(lib_B_sequences)
  
# Merging Data

gecko <- merge(gecko, lib_B_sequences, c("UID", "Gene"))

nrow(gecko)



# Computing Log2 Fold Change on BC3 data ----------------------------------

# Extracting raw counts

raw_counts <- as.matrix(gecko[, c("BC-3_Day0_Rep1", "BC-3_Day0_Rep2", 
                                      "BC-3_Day0_Rep3","BC-3_Day14_Rep1", 
                                      "BC-3_Day14_Rep2", "BC-3_Day14_Rep3")])

rownames(raw_counts) <- gecko$UID

condition_data <- data.frame(
  condition = factor(
    c("Day0", "Day0", "Day0","Day14", "Day14", "Day14"),
    levels = c("Day0", "Day14")
  )
)

rownames(condition_data) <- colnames(raw_counts)

# Running DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = condition_data,
  design = ~condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Day14", "Day0"))

# Extracting Log2 Fold Change data

gecko$LFC <- res[gecko$UID, "log2FoldChange"]

# Filter out infinite LFC values
valid <- is.finite(gecko$LFC)

gecko_model <- gecko[valid, , drop=FALSE]

# Setting dependent variable y
y <- gecko_model$LFC

# Visualizing Log2 Fold Change distribution
hist(y, breaks = 50, main = "LFC Distribution", xlab = "LFC")


# Computing Independent Variables -----------------------------------------

# Computing GC content

gecko_model$GC <- (str_count(gecko_model$Guide, "G") + 
                     str_count(gecko_model$Guide, "C")
                   )/str_length(gecko_model$Guide)

# Visualizing GC-content distribution

hist(gecko_model$GC, breaks = 20, main = "GC-content", xlab = "GC")

# Transforming sequence to numerical representations for each position

sequence_matrix <- str_split_fixed(
  gecko_model$Guide, 
  pattern = "", 
  n = 20)

X <- do.call(
  cbind,
  lapply(1:20, function(i) {
    model.matrix(~ sequence_matrix[, i] - 1)
  })
)

colnames(X) <- paste0(
  rep(paste0("pos", 1:20), each = 4),
  "_",
  rep(c("A", "C", "G", "T"), 20)
)

# Merging GC-content and position data

X <- cbind(X, GC = gecko_model$GC)




# Splitting Data ----------------------------------------------------------

set.seed(2026)

genes <- unique(gecko_model$Gene)


train_genes <- sample(genes, floor(0.8 * length(genes)))
test_genes <- setdiff(genes, train_genes)

train_index <- which(gecko_model$Gene %in% train_genes)
test_index <- which(gecko_model$Gene %in% test_genes)

# Creating DMatrix representations for XGBoost

dtrain <- xgb.DMatrix(
  X[train_index, , drop = FALSE],
  label = y[train_index]
)

dtest <- xgb.DMatrix(
  X[test_index, , drop = FALSE],
  label = y[test_index]
)


# XGBoost -----------------------------------------------------------------

model <- xgb.train(
  data = dtrain,
  nrounds = 1000,
  params = xgb.params(
    objective = "reg:squarederror",
    nthread = 1,
    max_depth = 4,
    eta = 0.02,
    min_child_weight = 10,
    subsample = 0.7,
    colsample_bytree = 0.8
  )
)

# Predicting with model

pred_train <- predict(model, dtrain)
pred_test <- predict(model, dtest)


cor(y[train_index], pred_train, method = "spearman")
cor(y[test_index], pred_test, method = "spearman")

#xgb.importance(model)
