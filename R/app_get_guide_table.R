# app_get_guide_table.R
# Code for computing extracting guide sequence and related metrics

# Purpose:
# Extracting sgRNAs and related metrics for a valid gene ID

# Input:
# - Sources R/app_finding_sgRNAs.R

# Output:
# - Guide data frame with guide, ID, strand and chromosome and attributes
# for gene ID and name

# Reading in Necessary Libraries ------------------------------------------

source("R/app_finding_sgRNAs.R")

# Extracting sgRNAs and Data ----------------------------------------------

get_guide_table <- function(geneID) {
  # Gene as GRange object
  gene <- genes(
    ensembl,
    filter = ~gene_id == geneID
  )
  
  # Changing style to match datasets
  seqlevelsStyle(gene) <- "UCSC"
  # Extracting GOI from reference genome --> DNAString Class
  gene_sequence <- getSeq(genome, gene)[[1]]
  # Creating reverese complement of GOI
  gene_sequence_reverse <- reverseComplement(gene_sequence)
  
  # Finding sgRNAs for GOI --------------------------------------------------
  forward_guides <- finding_sgRNAs(gene_sequence)
  forward_guides$strand <- as.character(strand(gene))
  
  reverse_guides <- finding_sgRNAs(gene_sequence_reverse)
  reverse_guides$strand <- ifelse(
    as.character(strand(gene)) == "+", "-", "+")
  
  # Creating data frame with all forward and reverse sgRNAs
  guides <- rbind(forward_guides, reverse_guides)
  
  # Adding ID to guides -----------------------------------------------------
  guides$ID <- 1:nrow(guides)
  
  # Adding Gene Name
  attr(guides,"gene_symbol") <- gene$gene_name[1]
  attr(guides, "ensembl_id") <- geneID
  
  # Adding Chromosome Information -------------------------------------------
  guides$chromosome <- as.character(seqnames(gene))
  
  # Constructing Genomic Ranges for Guides ----------------------------------
  gene_start <- start(gene)
  gene_end   <- end(gene)
  
  guides$genomic_start <- if (as.character(strand(gene)) == "+") {
    gene_start + guides$start - 1 
  } else { 
    gene_end - guides$end + 1 
  }
  
  guides$genomic_end <- if (as.character(strand(gene)) == "+") {
    gene_start + guides$end - 1
  } else {
    gene_end - guides$start + 1
  }
  
  return( guides )
}
