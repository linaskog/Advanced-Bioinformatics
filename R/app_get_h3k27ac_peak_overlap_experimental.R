# app_get_h3k27ac_peak_overlap.R

# Purpose:
# - Find h3K27ac peak overlaps for data frame of guides

# Input:
# - Sources R/app_get_gene_overlap_count.R
# - h3k27ac_peaks

# Output:
# - Guide data frame from app_get_gene_overlap_count.R with H3K27ac peak overlap
# column (TRUE/FALSE)

# Loading Libraries -------------------------------------------------------

h3k27ac_peaks <- rtracklayer::import(
  "~/Advanced-Bioinformatics/GSM4040945_BC3.merge.peak.bed.gz",
  format = "BED"
)
#source("R/app_get_gene_overlap_count.R")
source("R/app_get_gene_overlap_count_experimental.R")


# H3K27ac -----------------------------------------------------------------

get_h3k27ac_peak_overlap <- function(geneID) {

  # Reading in guid table from get_gene_overlap_count
  guides <- get_gene_overlap_count(geneID)

  # Create grances for guides
  guide_gr <- GRanges(
    seqnames = guides$chromosome,
    ranges = IRanges(
      start = guides$genomic_start,
      end = guides$genomic_end
    ),
    strand = guides$strand,
    UID = guides$ID)

  # Finding overlaps between guides and h3k27ac peaks
  guide_h3k27ac_hits <- findOverlaps(
    guide_gr,
    h3k27ac_peaks,
    ignore.strand = TRUE)

  # Adding colum for h3K27ac peak overlap (T/F)
  guides$H3K27ac_peak_overlap <- FALSE
  guides$H3K27ac_peak_overlap[unique(queryHits(guide_h3k27ac_hits))] <- TRUE

  return( guides )
}
