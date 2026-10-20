
# Reading in Datasets -----------------------------------------------------

# Human reference genome hg19
library(BSgenome.Hsapiens.UCSC.hg19)

# Human genome annotated with Ensembl IDs
library(EnsDb.Hsapiens.v75)


# Reading in necessary libraries ------------------------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)
library(xgboost)
library(stringr)


# Reading in model --------------------------------------------------------

xgb_model <- xgb.load("xgb_model.ubj")

# Finding genes -----------------------------------------------------------

INPUT <- "ENSG00000107485"

ensembl <- EnsDb.Hsapiens.v75

# Creating GRange class for gene of interest
gene <- genes(
  ensembl,
  filter = ~gene_id == INPUT
)

# Changing style to match datasets
seqlevelsStyle(gene) <- "UCSC"

# Extracting GOI from reference genome --> DNAString Class
gene_sequence <- getSeq(Hsapiens, gene)[[1]]

# Creating reverese complement of GOI
gene_sequence_reverse <- reverseComplement(gene_sequence)



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
  data.frame( "Guide" = as.character(extractAt(
               gene_sequence,
               guide_ranges
               )
               )
              )
}

# Creating data frame with all forward and reverse sgRNAs
guides <- rbind(finding_sgRNAs(gene_sequence), 
  finding_sgRNAs(gene_sequence_reverse))

# Adding identification column to guide data frame
guides$ID <- 1:nrow(guides)



# Guide Preprocessing -----------------------------------------------------

guides$GC <- (str_count(guides$Guide, "G") + 
                str_count(guides$Guide, "C")
)/str_length(guides$Guide)

sequence_matrix <- str_split_fixed(
  guides$Guide, 
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

X <- cbind(X, GC = guides$GC)


predictions <- predict(xgb_model, X)

