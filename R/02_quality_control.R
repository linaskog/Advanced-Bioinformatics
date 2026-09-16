# 02_quality_control.R
# Advanced Bioinformatics: sgRNA quality project
# Purpose:
# Assess the quality and consistency of the imported CRISPR screen data.
# Input:
# Objects created by R/01_import_data.R:
# - gecko_counts
# - gecko_library
# - bc3_mageck
# - bjab_mageck
#
# Checks performed:
# 1. Dataset dimensions and data types
# 2. Missing values
# 3. Duplicate and unmatched sgRNA identifiers
# 4. Zero and negative counts
# 5. Sequencing depth across samples
# 6. Agreement between biological replicates
#
# Output:
# Quality-control summaries and diagnostic plots.
# This script identifies potential data-quality problems.
# Import required data source("R/01_import_data.R")

source("R/01_import_data.R")

# check dataset dimensions
dataset_dimensions <- data.frame(
  dataset=c("gecko_counts",
            "gecko_library",
            "bc3_mageck", 
            "bjab_mageck")
  ,rows= c(nrow(gecko_counts),
           nrow(gecko_library),
           nrow(bc3_mageck),
           nrow(bjab_mageck))
           
  ,columns= c(ncol(gecko_counts),
              ncol(gecko_library),
              ncol(bc3_mageck),
              ncol(bjab_mageck))
)

dataset_dimensions

# the data type in the gecko_count file
count_column_type <- sapply(gecko_counts, class)
count_column_type

library_column_type<- sapply(gecko_library, class)
library_column_type

# Check data types in the MAGeCK tables 
bc3_mageck_types <- sapply(bc3_mageck, class)
bjab_mageck_types <- sapply(bjab_mageck, class)

bc3_mageck_types
bjab_mageck_types

identical(bc3_mageck_types, bjab_mageck_types)

# missing value
missing_count<- colSums(is.na(gecko_counts))
missing_count
sum(missing_count)

missing_summary<- data.frame(
  datasets=c(
    "gecko_counts", 
    "gecko_library",
    "bc3_mageck",
    "bjab_mageck"
    ),
  missing_values=c(
    sum(is.na(gecko_counts)),
    sum(is.na(gecko_library)),
    sum(is.na(bc3_mageck)),
    sum(is.na(bjab_mageck))
  )
)
missing_summary

# check for duplicate
duplicate_summary<- data.frame(
  dataset=c(
    "gecko_counts", 
    "gecko_library",
    "bc3_mageck",
    "bjab_mageck"
  ),
  identifier=c(
    "UID",
    "UID",
    "id",
    "id"
  ),
  duplicate_identifiers=c(
    sum(duplicated(gecko_counts$UID)),
    sum(duplicated(gecko_library$UID)),
    sum(duplicated(bc3_mageck$id)),
    sum(duplicated(bjab_mageck$id))
  )
)
duplicate_summary

# Check whether sgRNA identifiers match between tables
uid_missing_from_library<- setdiff(
  gecko_counts$UID,
  gecko_library$UID
  )
uid_missing_from_counts<- setdiff(
  gecko_library$UID,
  gecko_counts$UID
)
uid_match_summary <-data.frame(
  check= c(
    "Count UIDs missing from library",
    "Library UIDs missing from counts"
  ),
  number_missing= c(
    length(uid_missing_from_library),
    length(uid_missing_from_counts)
  )
)
uid_match_summary


# Check gene assignments for matching UIDs

library_gene_matched<- gecko_library$Gene[match(
  gecko_counts$UID, gecko_library$UID)
  ]
gene_assignment_mismatches<- sum(
  gecko_counts$Gene != library_gene_matched)
gene_assignment_mismatches


# Inspect examples of mismatched gene assignments 

mismatch_rows <- which(
  gecko_counts$Gene != library_gene_matched
)

gene_mismatch_examples <- data.frame(
  UID = gecko_counts$UID[mismatch_rows],
  gene_in_counts = gecko_counts$Gene[mismatch_rows],
  gene_in_library = library_gene_matched[mismatch_rows]
)

head(gene_mismatch_examples, 10)

#inspect the mismatch data
gecko_counts[
  gecko_counts$Gene == "hsa-let-7a-1",
  c("UID", "Gene")
]

gecko_library[
  grepl("let-7a", gecko_library$Gene, ignore.case = TRUE),
  c("UID", "Gene", "Guide")
]

# The count and sequence tables contain the same UID set, but 42,236 UIDs
# have different gene annotations. Some differences result from Excel date
# conversion, while others suggest inconsistent UID-to-gene mappings.

# Summarise matched and mismatched gene assignments
gene_assignment_status <- data.frame(
  UID= gecko_counts$UID,
  gene_in_counts= gecko_counts$Gene,
  gene_in_library= library_gene_matched,
  assignment_matches= gecko_counts$Gene == library_gene_matched
  )

assignment_summary<- data.frame(
  status=c("matched", "mismatched"),
  number_of_guides=c(
    sum(gene_assignment_status$assignment_matches),
    sum(!gene_assignment_status$assignment_matches)
    )
)
assignment_summary$percentage<- round(
  100* assignment_summary$number_of_guides/ nrow(gene_assignment_status),2 
)

assignment_summary

# Count genes before and after removing mismatched assignments
gene_lose_summary<- data.frame(
  datasets=c("all guides", "matched guides only"),
    unique_genes= c(
      length(unique(gene_assignment_status$gene_in_counts)),
      length(unique(gene_assignment_status$gene_in_counts[
        gene_assignment_status$assignment_matches]
    ) )
  )
)

gene_lose_summary

#So 1,867 targets can disappear, which is about 8.52% of all targets.
#the paper reports 1,864 miRNAs, so this similar number
# suggests that most completely lost targets may be miRNAs. 

# Identify targets that would disappear completely
all_genes <- unique(gene_assignment_status$gene_in_counts)

# Unique targets that have at least one correctly matched guide

matched_genes <- unique(
  gene_assignment_status$gene_in_counts[
  gene_assignment_status$assignment_matches]
  )

# Targets with no correctly matched guides
lost_genes<- setdiff(all_genes, matched_genes) 
length(lost_genes)
head(lost_genes,20)

# Identify miRNA targets among the completely lost targets
lost_mirna <- grepl("^hsa-(mir|let-)",
                    lost_genes,
                    ignore.case = TRUE)

lost_target_summary <- data.frame(
  category=c("miRNA targets","other targets"),
  numbers= c(
    sum(lost_mirna),
    sum(!lost_mirna)
  )
)
lost_target_summary

# Display any lost targets that are not recognised as miRNAs
lost_genes[!lost_mirna]

# All 1,864 miRNA targets would be lost if only matching UID-to-gene
# assignments were retained. Three additional targets (1-Dec, PLN and SPHAR)
# would also be lost. Therefore, the mismatches are systematic and primarily
# concern the miRNA section of the GeCKO v2 library, rather than random rows.

gecko_counts[
  gecko_counts$Gene %in% c("PLN", "SPHAR"),
  c("Gene", "UID")
]

gecko_library[
  gecko_library$Gene %in% c("PLN","SPHAR"),
  c("Gene","UID","Guide")
]

# Inspect the reciprocal UID mappings for PLN and SPHAR

uids_to_check <- c(
  "HGLibA_44493",
  "HGLibA_37205",
  "HGLibA_53929",
  "HGLibA_46641"
)

gecko_counts[
  gecko_counts$UID %in% uids_to_check,
  c("UID","Gene")
]

gecko_library[
  gecko_library$UID %in% uids_to_check,
  c("UID", "Gene", "Guide")
]


is_mirna_guides <- grepl(
  "^hsa-(mir|let-)",
  gene_assignment_status$gene_in_counts,
  ignore.case = TRUE
)

mismatched_guide_summary <- data.frame(
  category = c("miRNA guides", "Other guides"),
  mismatched_guides = c(
    sum(!gene_assignment_status$assignment_matches & is_mirna_guides),
    sum(!gene_assignment_status$assignment_matches & !is_mirna_guides)
  )
)

mismatched_guide_summary




# Keep non-miRNA targets
non_mirna_status <- gene_assignment_status[!is_mirna_guides, ]

# Count matched guides remaining for each target
matched_guides_per_gene <- aggregate(
  assignment_matches ~ gene_in_counts,
  data = non_mirna_status,
  FUN = sum
)

names(matched_guides_per_gene) <- c(
  "Gene",
  "number_of_matched_guides"
)

# Show the distribution
matched_guide_distribution <- as.data.frame(
  table(matched_guides_per_gene$number_of_matched_guides)
)

names(matched_guide_distribution) <- c(
  "number_of_matched_guides",
  "number_of_genes"
)

matched_guide_distribution
#7343 genes retain 6 matche guides
#11,600 genes retain exactly three matched guides  
#mapping problem may differ systematically between the A and B libraries.
# Determine whether each guide belongs to library A or library B

guide_library <- ifelse(
  grepl(
    "^HGLibA_", gene_assignment_status$UID
    ),
  "libraryA",
  ifelse(
    grepl("^HGLibB_", gene_assignment_status$UID
          ),
    "libraryB",
    "Others"
  )
  )

library_mapping_summary <- data.frame(
  library=c("libraryA","libraryB","Others"),
  total_guides=c(
    sum(guide_library== "libraryA"),
    sum(guide_library=="libraryB"),
    sum(guide_library=="Others")
  ),
  matched_guides=c(
    sum(gene_assignment_status$assignment_matches & guide_library== "libraryA"),
    sum(gene_assignment_status$assignment_matches & guide_library=="libraryB"),
    sum(gene_assignment_status$assignment_matches & guide_library=="Others")
  )
)

library_mapping_summary$mismatched_guides <- 
  library_mapping_summary$total_guides-
  library_mapping_summary$matched_guides

library_mapping_summary$matched_percentage <- round(
  100 * library_mapping_summary$matched_guides /
    library_mapping_summary$total_guides,
  2
)
library_mapping_summary

# Inspect the three mismatches in Library B

library_b_mismatches <- gene_assignment_status[
  guide_library == "libraryB" &
    !gene_assignment_status$assignment_matches,
]

library_b_mismatches

# library B contain 3 mismatch 1-DEC
# Create cleaned gene annotations without modifying the raw data
gene_counts_clean<- gecko_counts$Gene
gene_library_clean<- library_gene_matched

#Correct the Excel date conversion of DEC1
gene_counts_clean[gene_counts_clean== "1-Dec"] <- "DEC1"
gene_library_clean[gene_library_clean== "42339"] <- "DEC1"

clean_assignment_matches <-
  gene_counts_clean == gene_library_clean
# Confirm mismatches remaining in Library B
sum(
  guide_library == "libraryB" &
    !clean_assignment_matches
)

# Conclusion:
# Library B has no genuine UID-to-gene mismatches.
# Its three apparent mismatches were caused by Excel converting DEC1
# into "1-Dec" in one sheet and the date serial 42339 in the other.
# Library A still contains systematic mapping inconsistencies and should
# not be joined to guide sequences until the problem is resolved.


# Identify the experimental count columns
 count_columns<- setdiff(
   names(gecko_counts),
   c("UID","Gene")
 )

 # Count negative values in every sample
 negative_count_summary<- colSums(
   gecko_counts[count_columns]<0
 )
negative_count_summary 
sum(negative_count_summary)

# Count zero-read sgRNAs in each sample
zero_counts_summary<- data.frame(
  sample= count_columns,
  zero_counts_guides= colSums(
    gecko_counts[count_columns]==0
  )
)

zero_counts_summary$zero_count_precentage<-round(
  100*zero_counts_summary$zero_counts_guides/
    nrow(gecko_counts),
  2
)
zero_counts_summary

# 60% zero in BJAB replicate3 day 14
# Calculate the total sequencing depth of each sample

sequencing_depth_summary<- data.frame(
  sample=count_columns,
  total_reads= colSums(gecko_counts[count_columns])
)

sequencing_depth_summary$total_reads_millions<- round(
  sequencing_depth_summary$total_reads/1000000,
  2
)
sequencing_depth_summary

#comparing samples
# Log-transform counts before comparing samples
log2_counts<- log2(
  as.matrix(gecko_counts[count_columns])+1
)

# Calculate correlations between all samples
replicate_correlations<- cor(
  log2_counts,
  method = "pearson"
)
round(replicate_correlations, 2)


# Replicate-correlation interpretation:
# BC-3 Day 0 replicates show very strong agreement (r = 0.94-0.96).
# BC-3 Day 14 replicates show moderate-to-good agreement (r = 0.80-0.83).
# BJAB Day 0 replicates show strong agreement (r = 0.88-0.89).
# BJAB Day 14 replicates show weak agreement (r = 0.54-0.58).
# BJAB Day 14 also contains more zero-count guides. These samples should
# not be removed automatically, but their variability must be considered
# when constructing and evaluating the model.

