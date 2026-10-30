# app_get_gene_overlap_count.R
# code for finding overlap counts of sgRNA's

# Purpose:
# Find gene overlap counts for guides for a given gene ID

# Input:
# - Sources R/app_get_region_type.R
# - Bowtie index for hg19 reference genome, computed with SLURM

# Output:
# - Guide dataframe from app_get_region_type.R with gene overlap count column

# Reading in Necessary Libraries ------------------------------------------
library(Rbowtie)
#source("R/app_get_region_type.R")
source("R/app_get_region_type_experimental.R")
# Reading in Datasets -----------------------------------------------------
# Path to the hg19 Bowtie index
bowtie_index <- "hg19_index/hg19"
# Path for Bowtie alignment results


# Computing gene overlap count --------------------------------------------
get_gene_overlap_count <- function(geneID) {
  guides <- get_region_type(geneID)

  # Create temporary files for this gene analysis
  guide_fasta <- tempfile(fileext = ".fa")
  max_fasta <- tempfile(fileext = ".fa")
  alignment_file <- tempfile(fileext = ".tsv")

  # Remove temporary files when the function finishes or fails
  on.exit(
    unlink(c(guide_fasta, max_fasta, alignment_file)),
    add = TRUE
  )
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
    m=100,
    f = TRUE,
    max=max_fasta,
    outfile = alignment_file,
    force = TRUE
  )

  # Identify guides suppressed because they have more than 100 matches
  suppressed_ids <- character(0)

  if (file.exists(max_fasta)) {
    max_lines <- readLines(max_fasta, warn = FALSE)

    suppressed_ids <- sub(
      "^>",
      "",
      grep("^>", max_lines, value = TRUE)
    )
  }

  # Flag highly multi-mapping guides
  guides$high_multimapping <- paste0("guide_", guides$ID) %in% suppressed_ids


  # Read Bowtie alignments, including the case of no reported matches
  if (!file.exists(alignment_file) ||
      file.info(alignment_file)$size == 0) {

    # Create an empty table with the expected columns
    guide_coords <- data.frame(
      UID = character(),
      strand = character(),
      chromosome = character(),
      start_0 = integer(),
      guide = character(),
      quality = character(),
      mismatch = integer()
    )

  } else {

    guide_coords <- read.delim(
      alignment_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )

    if (ncol(guide_coords) != 8L) {
      stop("Unexpected Bowtie output format: expected 8 columns")
    }

    guide_coords <- guide_coords[, 1:7, drop = FALSE]

    colnames(guide_coords) <- c(
      "UID", "strand", "chromosome", "start_0",
      "guide", "quality", "mismatch"
    )
  }



  global_guide_gr <- GRanges(
    seqnames = guide_coords$chromosome,
    ranges = IRanges(
      start = guide_coords$start_0 + 1L,
      #this calculates the correct ending position from the actual aligned sequence length
      end = guide_coords$start_0 + nchar(guide_coords$guide)
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

  # Count distinct overlapping genes for each guide
  genes_per_guide <- table(guide_gene_pairs_unique$UID)

  # Match counts using guide IDs, not table positions
  overlap_counts <- as.integer(
    genes_per_guide[as.character(guides$ID)]
  )

  # Guides with no identified gene overlaps receive zero
  overlap_counts[is.na(overlap_counts)] <- 0L

  # Highly multi-mapping guides have unknown exact overlap counts
  overlap_counts[guides$high_multimapping] <- NA_integer_

  # Add the corrected counts to the guide table
  guides$gene_overlap_count <- overlap_counts

  return( guides )
}
