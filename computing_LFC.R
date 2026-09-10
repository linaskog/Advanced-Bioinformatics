
# Reading in Packages -----------------------------------------------------

library(readxl)
library(tidyverse)
library(stringr)
library(DESeq2)


# Reading in Data ---------------------------------------------------------

excel_sheets("41467_2018_5506_MOESM4_ESM.xlsx")

# Raw Count Data
gecko <- read_excel("41467_2018_5506_MOESM4_ESM.xlsx", sheet = "C. GECKO_BC3_BJAB" )

# sgRNA Sequences
sequences <- read_excel("41467_2018_5506_MOESM4_ESM.xlsx", sheet = "F. GeCKOlibrary2")


# Merging Data
new_gecko <- merge(gecko, sequences, c("UID"))



# Computing Log2 Fold Change ----------------------------------------------

raw_counts <- as.matrix(new_gecko[, c("BC-3_Day0_Rep1", "BC-3_Day0_Rep2", 
                                      "BC-3_Day0_Rep3","BC-3_Day14_Rep1", 
                                      "BC-3_Day14_Rep2", "BC-3_Day14_Rep3")])
rownames(raw_counts) <- new_gecko$UID

condition_data <- data.frame(
  condition = factor(
    c("Day0", "Day0", "Day0","Day14", "Day14", "Day14"),
    levels = c("Day0", "Day14")
  )
)

rownames(condition_data) <- colnames(raw_counts)

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = condition_data,
  design = ~condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Day14", "Day0"))

new_gecko$LFC <- res[new_gecko$UID, "log2FoldChange"]


