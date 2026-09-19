# 03_data_preprocessing.R
# Advanced Bioinformatics: sgRNA quality project

# Purpose:
# Process data prior to model training
# Independent variable: log2 fold change of raw sgRNA counts for BC3 data
# Dependent variables: GC-content and position base data for sgRNA guides

# Reading in Packages -----------------------------------------------------
library(readxl)
library(stringr)
library(DESeq2)


# Reading in Data ---------------------------------------------------------
filepath = "~/Advanced-Bioinformatics/41467_2018_5506_MOESM4_ESM.xlsx"
excel_sheets(filepath)

# Raw Count Data 
gecko_data <- read_excel(filepath, sheet = "C. GECKO_BC3_BJAB" )

# Library Information
gecko_library <- read_excel(filepath, sheet = "F. GeCKOlibrary2")

# Extracting Library B 
libB_index <- grep("LibB", gecko_library$UID)
libB <- gecko_library[libB_index, ]

print(paste("Number of sgRNAs in library B:", nrow(libB))) #(58028 sgRNAs)

# Removing Non-gene-targeting sgRNAs
controls <- grep("NonTargeting", libB$Gene)
clean_libB <- libB[-controls, ]

print(paste("Number of gene-targeting sgRNAs in Library B:", 
            nrow(clean_libB))) # (57028 sgRNAs)

# Merging Data
gecko <- merge(gecko_data, clean_libB, c("UID", "Gene"))

print(paste("Number of remaining sgRNAs:", nrow(gecko))) #(57025 sgRNAs)



# Computing Log2 Fold Change on BC3 data ----------------------------------

# Extracting raw counts of BC3 data
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

# Setting dependent variable y
y <- gecko$LFC

# Visualizing Log2 Fold Change distribution
hist(y, breaks = 50, main = "LFC Distribution", xlab = "LFC")


# Computing Independent Variables -----------------------------------------

# Computing GC content

gecko$GC <- (str_count(gecko$Guide, "G") + 
               str_count(gecko$Guide, "C")
)/str_length(gecko$Guide)

# Visualizing GC-content distribution

hist(gecko$GC, breaks = 20, main = "GC-content", xlab = "GC")

# Transforming sequence to numerical representations for each position

sequence_matrix <- str_split_fixed(
  gecko$Guide, 
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


# Splitting Data ----------------------------------------------------------

set.seed(2026)

genes <- unique(gecko$Gene)


train_genes <- sample(genes, floor(0.8 * length(genes)))
test_genes <- setdiff(genes, train_genes)

train_index <- which(gecko$Gene %in% train_genes)
test_index <- which(gecko$Gene %in% test_genes)

# Creating DMatrix representations for XGBoost

dtrain <- xgb.DMatrix(
  X[train_index, , drop = FALSE],
  label = y[train_index]
)

dtest <- xgb.DMatrix(
  X[test_index, , drop = FALSE],
  label = y[test_index]
)
