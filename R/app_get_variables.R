# app_get_variables.R

# Purpose:
# Prepare variable data for modelling

# Input:
# - sources app

# Sourcing extracted data -------------------------------------------------

source("~/GitHub/Advanced-Bioinformatics/R/app_get_h3k27ac_peak_overlap.R")

# Preparing Independent Variables -----------------------------------------
get_variables <- function(geneID) {
  guides <- get_h3k27ac_peak_overlap(geneID)
  
  variables <- data.frame("UID" = guides$ID,
                          "Guide" = guides$guide)
  
  # Extracting Position-Wise Nucleotide Information
  sequence_matrix <- str_split_fixed(
    guides$guide, 
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
                     GC = (str_count(guides$guide, "G") + 
                             str_count(guides$guide, "C")) / 
                       str_length(guides$guide)
  )
  
  # Adding Strand Data to Variables
  variables <- cbind(variables, strand = as.numeric(factor(guides$strand,
                                                           levels = c("-", "+"))))
  
  # Adding Chromosome Information to Variables
  variables <- cbind(variables, 
                     chromosome = as.numeric(
                       factor(
                         guides$chromosome, 
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
                         guides$region_type, 
                         levels = c("exon", "intergenic", "intron")
                       )
                     )
  )
  
  # Adding Gene Overlap Count
  variables <- cbind(variables, 
                     gene_overlap_count = guides$gene_overlap_count)
  
  # Adding H3K27ac Peak Overlap (True = 1, False = 0)
  variables <- cbind(variables, 
                     H3K27ac_peak_overlap = as.numeric(
                       guides$H3K27ac_peak_overlap)
  )
}
