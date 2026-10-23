# app_get_region_type.R
# Function for finding region type of sgRNAs

# Purpose:
# Find region type associated with guide sequence
# - exon, intron or intergenic

# Input:
# - Sources R/app_get_guide_table.R
# - hg19 reference genome as gtf

# Output:
# - guide data frame from app_get_guide_table.R with added region type column

# Reading in Necessary Libraries ------------------------------------------

source("~/GitHub/Advanced-Bioinformatics/R/app_get_guide_table.R")

# Reading in Datasets -----------------------------------------------------

gtf <- rtracklayer::import(
  "~/GitHub/Advanced-Bioinformatics/data/reference/hg19.refGene.gtf.gz"
)

transcripts<- gtf[gtf$type== "transcript"]
exons <- gtf[gtf$type == "exon"]

get_region_type <- function(geneID) {
  
  guides <- get_guide_table(geneID)
  
  guide_gr <- GRanges(
    seqnames = guides$chromosome,
    ranges = IRanges(
      start = guides$genomic_start,
      end = guides$genomic_end
    ),
    strand = guides$strand,
    UID = guides$ID)
  
  # Adding Region Type Information ------------------------------------------
  
  guide_transcript_hits <- findOverlaps(
    guide_gr,
    transcripts,
    ignore.strand = TRUE
  )
  
  # guides overlapping exons
  guide_exon_hits <- findOverlaps(
    guide_gr,
    exons,
    ignore.strand = TRUE
  )
  
  # guides overlapping introns
  transcript_guides <- unique(queryHits(guide_transcript_hits))
  exon_guides <- unique(queryHits(guide_exon_hits))
  intron_only_guides <- setdiff(transcript_guides, exon_guides)
  
  guides$region_type <- "intergenic"
  guides$region_type[intron_only_guides] <- "intron"
  guides$region_type[exon_guides] <- "exon"
  
  return( guides )
}
