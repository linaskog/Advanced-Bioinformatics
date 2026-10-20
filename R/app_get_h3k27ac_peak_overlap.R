# app_get_h3k27ac_peak_overlap.R

# Purpose:

# Input:


# Loading Libraries -------------------------------------------------------

h3k27ac_peaks <- rtracklayer::import(
  "~/Advanced-Bioinformatics/GSM4040945_BC3.merge.peak.bed.gz",
  format = "BED"
)
source("~/GitHub/Advanced-Bioinformatics/R/app_get_gene_overlap_count.R")

# H3K27ac -----------------------------------------------------------------

get_h3k27ac_peak_overlap <- function(geneID) {
  
  guides <- get_gene_overlap_count(geneID)
  
  guide_gr <- GRanges(
    seqnames = guides$chromosome,
    ranges = IRanges(
      start = guides$genomic_start,
      end = guides$genomic_end
    ),
    strand = guides$strand,
    UID = guides$ID)
  
  guide_h3k27ac_hits <- findOverlaps(
    guide_gr,
    h3k27ac_peaks,
    ignore.strand = TRUE)
    
  guides$H3K27ac_peak_overlap <- FALSE
    
  guides$H3K27ac_peak_overlap[unique(queryHits(guide_h3k27ac_hits))] <- TRUE
  
  return( guides )
}
  