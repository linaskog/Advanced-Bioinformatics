# 01_import_data.R
# Advanced Bioinformatics - sgRNA efficiency project
# Purpose:
# Import and inspect the CRISPR screen data from Manzano et al. (2018).
# Initial focus: GeCKO v2 screen in the BJAB cell line.
install.packages("readxl")
library(readxl)
library(tidyverse)
# Path to the supplementary CRISPR screen data
screen_file <- "data/raw/41467_2018_5506_MOESM4_ESM.xlsx"
excel_sheets(screen_file)

# import raw sgRNA counts for GecKo screen
gecko_counts<- read_excel(screen_file, sheet = "C. GECKO_BC3_BJAB")
dim(gecko_counts)
names(gecko_counts)

# Inspect sgRNAs targeting CTDP1 gene in both BC_3 nad BJAB
ctdp1_count<- gecko_counts[gecko_counts$Gene== "CTDP1",]
ctdp1_count

#calculate mean across the 3 BC3 replicates
ctdp1_count$BC3_Day0_mean<- rowMeans(ctdp1_count[, c(
  "BC-3_Day0_Rep1","BC-3_Day0_Rep2","BC-3_Day0_Rep3")]
)

ctdp1_count$BC3_Day14_mean<- rowMeans(ctdp1_count[, c(
  "BC-3_Day14_Rep1","BC-3_Day14_Rep2","BC-3_Day14_Rep3")]
)
ctdp1_count$BC3_Day0_mean 
ctdp1_count$BC3_Day14_mean
ctdp1_count[,c("UID", "Gene", "BC3_Day0_mean", "BC3_Day14_mean")]

# calculate log2fold change for CTDP1
ctdp1_count$BC3_log2FC<- log2(
  ctdp1_count$BC3_Day14_mean/ctdp1_count$BC3_Day0_mean
)
ctdp1_count[, c("UID","Gene","BC3_Day0_mean", "BC3_Day14_mean","BC3_log2FC")]


# Import GeCKO v2 guide sequences
gecko_library<- read_excel(screen_file, sheet = "F. GeCKOlibrary2")
dim(gecko_library)
names(gecko_library)

gecko_library[gecko_library$UID=="HGLibB_11592",]

# Add guide sequence to the CTDP1 count table
ctdp1_count$Guide<- gecko_library$Guide[
  match(ctdp1_count$UID, gecko_library$UID)
]
ctdp1_count[, c(
  "UID","Gene","BC3_Day0_mean", "BC3_Day14_mean","BC3_log2FC","Guide"
  )]


# Calculate GC content for each guide seq
# Calculate GC content for each guide sequence
ctdp1_count$GC_content <- sapply(
  ctdp1_count$Guide,
  function(seq) {
    bases <- strsplit(seq, "")[[1]]
    sum(bases %in% c("G", "C")) / length(bases)
  }
)

ctdp1_count[, c( "UID",
                 "Guide",
                 "GC_content",
                 "BC3_log2FC")]
