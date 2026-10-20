# app_get_gene_overlap_count.R

# Purpose:
# Find gene overlap counts

# Input:
# - Sources R/app_get_region_type.R

# Reading in Necessary Libraries ------------------------------------------
library(Rbowtie)
source("~/GitHub/Advanced-Bioinformatics/R/app_get_region_type.R")

# Reading in Datasets -----------------------------------------------------
bowtie_index <- "~/GitHub/Advanced-Bioinformatics/hg19_index/hg19"
alignment_file <- "~/GitHub/Advanced-Bioinformatics/data/processed/new_guides_hg19.tsv"

# Computing gene overlap count --------------------------------------------
get_gene_overlap_count <- function(geneID) {
  guides <- get_region_type(geneID)
  guide_fasta <- tempfile(fileext = ".fa")
  writeLines(
    unlist(Map(
      function(id, seq) c(paste0(">guide_", id), seq),
      guides$ID, guides$guide)),
    guide_fasta
  )
  
  Rbowtie::bowtie(
    sequences = guide_fasta,
    index = path.expand(bowtie_index),
    v = 0,
    a = TRUE,
    f = TRUE,
    outfile = alignment_file,
    force = TRUE
  )
  
  # Read alignments
  guide_coords <- read.delim(
    alignment_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(guide_coords) <- c("UID", "strand", "chromosome", "start_0",
                              "guide", "quality", "mismatch")
  
  global_guide_gr <- GRanges(
    seqnames = guide_coords$chromosome,
    ranges = IRanges(
      start = guide_coords$start_0 + 1L,
      end = guide_coords$start_0 + 21L 
    ),
    strand = guide_coords$strand,
    UID= as.integer(sub("^guide_", "", guide_coords$UID))
  )
  
  global_transcript_hits <- findOverlaps(
    global_guide_gr,
    transcripts,
    ignore.strand=TRUE
  )
  
  # Which gene does each guide overlap
  guide_gene_pairs <- data.frame(
    UID = global_guide_gr$UID[queryHits(global_transcript_hits)],
    gene_name= transcripts$gene_name[subjectHits(global_transcript_hits)],
    stringsAsFactors = FALSE
  )
  
  guide_gene_pairs_unique <- unique(guide_gene_pairs)
  
  genes_per_guide<- table(guide_gene_pairs_unique$UID)
  
  guides$gene_overlap_count <- as.integer(
    genes_per_guide[guides$ID])
  
  return( guides )
}
