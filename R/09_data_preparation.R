# 09_data_preparation
# Advanced Bioinformatics - sgRNA efficiency project
# Purpose: Prepare CRISPR screen data from QC, genome annotation and DeepSeq
#          data for modelling.
# Input: 
# - CRISPR screen data from GeCKO library b
# - Genome annotation data
# - ChIP-seq peak data annotations
#
# Output: Variable dataframe as RDS object containing:
# - guide ID (UID)
# - Gene name
# - Positionwise Nucleotide Data (80 columns)
# - GC-Ratio in guide
# - Strand information
# - Region type infromation (levels exon, intergenic, intron)
# - Chromosome information (levels chr1-chr22, chrX, chrY)
# - Gene overlap counts
# - H3K27ac peak overlap data

# Loading Necessary Libraries ---------------------------------------------

library(DESeq2)
library(ggplot2)
library(stringr)
library(tidyverse)

# Reading in Data ---------------------------------------------------------

library_b <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/library_b_prepared.rds")
library_b <- library_b[library_b$target_type == "targeting" & 
                                   library_b$sequence_status == "unique_sequence", ]
H3K27ac <- read_tsv(
  "~/GitHub/Advanced-Bioinformatics/data/processed/guide_hg19_bc3_h3k27ac.tsv")

# Merging Datasets
data <- merge(library_b, H3K27ac, by = "UID")

# Extracting variables of interest
data <- data[, c("UID", 
                 "Gene", 
                 "Guide", 
                 "BC-3_Day0_Rep1",
                 "BC-3_Day0_Rep2", 
                 "BC-3_Day0_Rep3",
                 "BC-3_Day14_Rep1", 
                 "BC-3_Day14_Rep2", 
                 "BC-3_Day14_Rep3",
                 "strand",
                 "chromosome",
                 "region_type",
                 "gene_overlap_count",
                 "H3K27ac_peak_overlap"
                 )
             ]


# Preparing Dependent Variable --------------------------------------------

# Preparing Dataframe for Storage of Variables
variables <- data.frame("UID" = data$UID,
                        "Gene" = data$Gene,
                        "Guide" = data$Guide)

# Extracting Raw sgRNA Count Data for BC-3 Cell-Line
raw_counts <- as.matrix(data[, grep("BC-3", colnames(data))])
rownames(raw_counts) <- data$UID

# Preparing Condition Data
condition_data <- data.frame(
  condition = factor(
    c("Day0", "Day0", "Day0", "Day14", "Day14", "Day14"),
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

res <- results(dds, 
               contrast = c("condition", "Day14", "Day0"))

# Saving Log2 Fold Change Variable
variables$LFC <- res[variables$UID, "log2FoldChange"]

# Computing Gene Average Log2 Fold Change
y_gene_avg <- tapply(variables$LFC, variables$Gene, mean, na.rm = TRUE)

# Centering Log2 Fold Change Variable 
variables$centered_LFC <- variables$LFC - y_gene_avg[variables$Gene]

# Visualizing Difference Between Raw LFC and Centered LFC 
ggplot(variables, aes(LFC)) + 
  geom_histogram(col = "black", fill = "lightgrey", binwidth = 0.1) +
  labs(
    title ="Raw LFC Histogram",
    x = "Log2 Fold Change"
  ) + 
  theme_bw()

ggplot(variables, aes(centered_LFC)) + 
  geom_histogram(col = "black", fill = "lightgrey", binwidth = 0.1) +
  labs(
    title ="Gene Centered LFC Histogram",
    x = "Centred Log2 Fold Change"
  ) + 
  theme_bw()

# Preparing Independent Variables -----------------------------------------

# Extracting Position-Wise Nucleotide Information
sequence_matrix <- str_split_fixed(
  data$Guide, 
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

variables <- cbind(variables, X)


# Adding Guide GC-Content to Variables
variables <- cbind(variables, 
                   GC = (str_count(data$Guide, "G") + 
                           str_count(data$Guide, "C")) / 
                     str_length(data$Guide)
)

# Adding Strand Data to Variables
variables <- cbind(variables, strand = as.numeric(factor(data$strand,
                                                         levels = c("-", "+"))))

# Adding Chromosome Information to Variables
variables <- cbind(variables, 
                   chromosome = as.numeric(
                     factor(
                       data$chromosome, 
                       levels = c( "chr1", "chr2", "chr3", "chr4", "chr5", 
                                   "chr6", "chr7", "chr8", "chr9", "chr10", 
                                   "chr11", "chr12", "chr13", "chr14", "chr15", 
                                   "chr16", "chr17", "chr18", "chr19", "chr20",
                                   "chr21", "chr22", "chrX", "chrY" )
                     )
                   )
)

# Adding region Type (levels "exon", "intergenic", "intron") to Variables

variables <- cbind(variables, 
                   region_type = as.numeric(
                     factor(
                       data$region_type, 
                       levels = c("exon",
                                  "intergenic",
                                  "intron")
                     )
                   )
                   )

# Adding Gene Overlap Count
variables <- cbind(variables, 
                   gene_overlap_count = data$gene_overlap_count)

# Adding H3K27ac Peak Overlap (True = 1, False = 0)
variables <- cbind(variables, 
                   H3K27ac_peak_overlap = as.numeric(
                     data$H3K27ac_peak_overlap)
                   )


# Saving Variables to RDS object ------------------------------------------
saveRDS(variables, "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")
