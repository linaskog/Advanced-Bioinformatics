# 01_import_data.R
# Advanced Bioinformatics: sgRNA quality project
# Purpose:
# Import the CRISPR screen data from Manzano et al. (2018).
# Input:
# 1. MOESM4: raw sgRNA counts and sgRNA library sequences
# 2. MOESM6: processed gene-level MAGeCK results
# Output:
# R objects containing the imported raw datasets.
#
# Important:
# This script only imports and verifies data.
# Load package library(readxl)
library(readxl)
screen_file <- "data/raw/41467_2018_5506_MOESM4_ESM.xlsx"
mageck_file <- "data/raw/41467_2018_5506_MOESM6_ESM.xlsx"

file.exists(screen_file)
file.exists(mageck_file)

screen_sheets<-excel_sheets(screen_file)
mageck_sheets <-excel_sheets(mageck_file)
screen_sheets
mageck_sheets

# Import raw GeCKO v2 screen counts 
gecko_counts <- read_excel(
  screen_file,
  sheet ="C. GECKO_BC3_BJAB" )
dim(gecko_counts)
names(gecko_counts)
head(gecko_counts)

# Import GeCKO v2 guide sequences
gecko_library<- read_excel(
  screen_file,
  sheet = "F. GeCKOlibrary2"  )
dim(gecko_library)
names(gecko_library)
head(gecko_library)

# Import gene-level MAGeCK results
mageck_sheets
bc3_mageck<- read_excel(mageck_file, sheet = "GECKO_BC-3")
bjab_mageck<-read_excel(mageck_file, sheet ="GECKO_BJAB")
dim(bc3_mageck)
dim(bjab_mageck)
names(bc3_mageck)
names(bjab_mageck)

# validate imported data
stopifnot(
  all(c("UID","Gene") %in% names(gecko_counts)),
  all(c("UID","Guide","Gene") %in% names(gecko_library)),
  all(c("id", "neg|fdr", "neg|lfc") %in% names(bc3_mageck)),
  all(c("id", "neg|fdr", "neg|lfc") %in% names(bjab_mageck))
)
message("all required datasets were imported successfully")

