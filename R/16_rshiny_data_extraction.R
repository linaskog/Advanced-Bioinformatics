# 16_rshiny_data_extraction.R
# Code for computing and extracting variables associated with guide sgRNAs

# Purpose:
# Finding sgRNA's and extracting necessary information for modelling


# Input:
# - Human reference genome
# - Annotated human genome
# - h3k27ac peak data
# - hg19 bowtie indexes


# Reading in Necessary Libraries ------------------------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)
library(stringr)
library(Rbowtie)

# Reading in Datasets -----------------------------------------------------

# Human reference genome hg19
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# Human genome annotated with Ensembl IDs
library(EnsDb.Hsapiens.v75)
ensembl <- EnsDb.Hsapiens.v75

gtf <- rtracklayer::import(
  "data/reference/hg19.refGene.gtf.gz"
)

h3k27ac_peaks <- rtracklayer::import(
  "~/Advanced-Bioinformatics/GSM4040945_BC3.merge.peak.bed.gz",
  format = "BED"
)

transcripts<- gtf[gtf$type== "transcript"]
exons <- gtf[gtf$type == "exon"]




# Finding sgRNAs ----------------------------------------------------------

finding_sgRNAs <- function(gene_sequence, PAM = "NGG") {
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



# Finding Example Gene ----------------------------------------------------

# Example input
INPUT <- "ENSG00000107485"

# Gene as GRange object
gene <- genes(
  ensembl,
  filter = ~gene_id == INPUT
)

# Changing style to match datasets
seqlevelsStyle(gene) <- "UCSC"
# Extracting GOI from reference genome --> DNAString Class
gene_sequence <- getSeq(genome, gene)[[1]]
# Creating reverese complement of GOI
gene_sequence_reverse <- reverseComplement(gene_sequence)

gene_start <- start(gene)
gene_end <- end(gene)


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

guide_gr <- GRanges(
  seqnames = as.character(seqnames(gene)),
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


# Adding Gene Overlap Count Information -----------------------------------

# Write sgRNAs to FASTA
guide_fasta <- tempfile(fileext = ".fa")
writeLines(
  unlist(Map(
    function(id, seq) c(paste0(">guide_", id), seq),
    guides$ID, guides$guide
    )
    ),
  guide_fasta
)
bowtie_index <- "~/GitHub/Advanced-Bioinformatics/hg19_index/hg19"
alignment_file <- "data/processed/new_guides_hg19.tsv"

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

#Which gene does each guide overlap
guide_gene_pairs <- data.frame(
  UID = guide_gr$UID[queryHits(global_transcript_hits)],
  gene_name= transcripts$gene_name[subjectHits(global_transcript_hits)],
  stringsAsFactors = FALSE
)

guide_gene_pairs_unique <- unique(guide_gene_pairs)

genes_per_guide<- table(guide_gene_pairs_unique$UID)

guides$gene_overlap_count <- as.integer(
  genes_per_guide[guides$ID])

# H3K27ac -----------------------------------------------------------------

guide_h3k27ac_hits <- findOverlaps(
  guide_gr,
  h3k27ac_peaks,
  ignore.strand = TRUE
)

guides$H3K27ac_peak_overlap <- FALSE

guides$H3K27ac_peak_overlap[
  unique(queryHits(guide_h3k27ac_hits))] <- TRUE
