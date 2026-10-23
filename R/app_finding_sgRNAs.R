# app_finding_sgRNAs.R
# Code for function extracting guide sequences from human reference genome 
# given a valid gene ID

# Purpose:
# Finding sgRNA's for input gene

# Input:
# - Ensembl geneID
# - Human reference genome
# - Ensembl nnotated human genome

# Output:
# - data frame with sgRNAs

# Reading in Necessary Libraries ------------------------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)

# Human reference genome hg19
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# Human genome annotated with Ensembl IDs
library(EnsDb.Hsapiens.v75)
ensembl <- EnsDb.Hsapiens.v75

# Finding sgRNAs ----------------------------------------------------------

finding_sgRNAs <- function(gene_sequence, 
                           PAM = "NGG") {
  # Locate PAM sequences in gene of interest
  pam_sequences <- matchPattern(
    PAM,
    gene_sequence,
    fixed = FALSE
  )
  
  # Locating PAM start positions for valid (20 nt) sgRNAs
  pam_start <- start(pam_sequences)
  pam_start <- pam_start[pam_start > 20]
  
  # Creating IRange object for each sgRNA
  guide_ranges <- IRanges(
    start = pam_start - 20,
    end = pam_start -1
  )
  
  # Creating data frame for Guide sequences extracted using IRanges
  data.frame(
    guide = as.character(extractAt(gene_sequence,
                                   guide_ranges)),
    start = start(guide_ranges),
    end = end(guide_ranges) )
}
