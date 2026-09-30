#06_sgRNA preprocessing


# Loading required packages -----------------------------------------------

library(stringr)
library(xgboost)

# Reading in data ---------------------------------------------------------

# This will be changed to tool that reads in sgRNA with IDs and guides
INPUT <- data.frame()

# Loading model
xgb_model <- xgb.load("xgb_model.ubj")


# Preparing data for modelling --------------------------------------------

sgRNA <- INPUT


# Extracing GC content
sgRNA$GC <- (str_count(sgRNA$Guide, "G") + str_count(sgRNA$Guide, "C")) / 
  str_length(sgRNA$Guide)


# One-hot encoding of sequences
sequence_matrix <- str_split_fixed(
  sgRNA$Guide, 
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
X <- cbind(X, GC = gecko$GC)



# Modelling ---------------------------------------------------------------

predictions <- predict(xgb_model, X)
