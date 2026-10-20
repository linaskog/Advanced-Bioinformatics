source("R/01_import_data.R")

# Identify whether each sgRNA belongs to Library A or Library B

guide_library<- ifelse(
  startsWith(gecko_counts$UID, "HGLibA_"),
  "libraryA",
  ifelse(
    startsWith(gecko_counts$UID, "HGLibB_"),
    "libraryB",
    "unknown"
  )
)

table(guide_library)     # A=65383    B=58028

#keep only the reliable Library B guides
library_b_counts<- gecko_counts[guide_library== "libraryB",]
dim(library_b_counts)
head(library_b_counts[,1:3])

# retain the library B guide seq 
library_b_sequences <- gecko_library[
  startsWith(gecko_library$UID, "HGLibB_"),
]
dim(library_b_sequences)
head(library_b_sequences)

# Join Library B counts and sequences by UID
# Rename the two Gene columns so we can compare them after joining
library_b_data <-dplyr::left_join(
  dplyr::rename(
    library_b_counts,
    Gene_counts= Gene
  ),
  dplyr::rename(
    library_b_sequences,
    Gene_library=Gene
  ),
  by= "UID"
)
dim(library_b_data)  
names(library_b_data)
head(library_b_data[, c("UID","Gene_counts","Gene_library","Guide")])

# confirm that every sgRNA got a seq
sum(is.na(library_b_data$Guide))

# Compare gene annotations before correcting DEC1
table(library_b_data$Gene_counts == 
        library_b_data$Gene_library,
      useNA = "ifany")

#library B annotation mismatch
library_b_mismatches <- library_b_data[
  library_b_data$Gene_counts !=
    library_b_data$Gene_library,
  c(
    "UID",
    "Gene_counts",
    "Gene_library",
    "Guide"
  )
]
library_b_mismatches

# Create cleaned gene annotations without modifying the raw data
library_b_data<- dplyr::mutate(
  library_b_data,
  Gene_counts_clean= ifelse(
    Gene_counts==  "1-Dec",
    "DEC1",
    Gene_counts
  ),
  Gene_library_clean= ifelse(
    Gene_library== "42339",
    "DEC1",
    Gene_library
  )
)

table(library_b_data$Gene_counts_clean== library_b_data$Gene_library_clean)

# Create the final gene annotation used in the analysis
library_b_data<- dplyr::mutate(
  library_b_data,
  Gene= Gene_counts_clean
)
sum(is.na(library_b_data$Gene))
length(unique(library_b_data$Gene))

# sgRNA sequence lengths
table(nchar(library_b_data$Guide))
# Check that sequences contain only A, C, G and T
table(grepl(
  "^[ACGT]+$",
  library_b_data$Guide))

# check for control and miRNA
mirna_candidates<- sort(
  unique(
    library_b_data$Gene[
    grepl(
      "mirna",
      library_b_data$Gene,
      ignore.case = TRUE
    )
  ]
 )
)
length(mirna_candidates)

control_candidates<- sort(
  unique(
    library_b_data$Gene[
      grepl(
        "control|non.?target",
        library_b_data$Gene,
        ignore.case = TRUE
      )
    ]
  )
)

length(control_candidates)
head(control_candidates,20)

#classification of targeting and non_targeting controls guides
library_b_data<- dplyr::mutate(
  library_b_data,
  target_type= ifelse(
    grepl(
      "^NonTargetingControlGuideForHuman_",
      Gene
    ),
    "non_targeting_control",
    "targeting"
  )
)

table(library_b_data$target_type)

# Create separate datasets without deleting anything

library_b_targeting<- library_b_data[library_b_data$target_type=="targeting",]
library_b_controls<- library_b_data[
  library_b_data$target_type=="non_targeting_control",]
dim(library_b_targeting)
dim(library_b_controls)
length(unique(library_b_targeting$Gene))

count_columns<- names(gecko_counts)[3:14]
sum(duplicated(library_b_data$UID))
sum(duplicated(library_b_data$Guide))
sum(is.na(library_b_data$Guide))
# Check for missing experimental counts
sum(is.na(library_b_data[count_columns]))

#1159 duplicated seq
# Collect every row whose guide sequence occurs more than once
duplicated_sequence_rows <- library_b_data |>
  dplyr::group_by(Guide) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::ungroup() |>
  dplyr::arrange(Guide)

# nr of rows involved in sequence duplication
nrow(duplicated_sequence_rows)

# nr of distinct repeated sequences
dplyr::n_distinct(duplicated_sequence_rows$Guide)
# Are repeated sequences targeting guides or controls?
table(duplicated_sequence_rows$target_type)

# Collect every row with the same guide seq
duplicated_sequence_rows<- library_b_data|> 
  dplyr::group_by(Guide)|>
  dplyr::filter(dplyr:: n() > 1) |>
  dplyr::ungroup()|>
  dplyr::arrange(Guide)
# nr of rows involved in sequence duplication
nrow(duplicated_sequence_rows)


# nr of distinct repeated sequences
dplyr::n_distinct(duplicated_sequence_rows$Guide)

#are they targeting guides or controls?
table(duplicated_sequence_rows$target_type)

# Summarize each repeated sequence
duplicated_sequence_summary <- duplicated_sequence_rows |>
  dplyr::group_by(Guide) |>
  dplyr::summarise(
    number_of_rows = dplyr::n(),
    number_of_genes = dplyr::n_distinct(Gene),
    .groups = "drop"
  )

# Check whether repeated sequences belong to one or multiple annotations
table(duplicated_sequence_summary$number_of_genes)

head(duplicated_sequence_rows)

#count-profile identifier for every repeated-sequence row
duplicated_sequence_rows$count_profile <- apply(
  as.data.frame(
    duplicated_sequence_rows[count_columns]
  ),
  1,
  paste,
  collapse= "|"
)

#different count profiles for each seque
duplicated_count_summary <- duplicated_sequence_rows |>
  dplyr::group_by(Guide) |>
  dplyr::summarise(
    number_of_rows = dplyr::n(),
    number_of_genes = dplyr::n_distinct(Gene),
    number_of_count_profiles =
      dplyr::n_distinct(count_profile),
    .groups = "drop"
  )
table(duplicated_count_summary$number_of_count_profiles)

example_sequence <- duplicated_count_summary$Guide[
  which.max(duplicated_count_summary$number_of_genes)
]

duplicated_sequence_rows |>
  dplyr::filter(Guide == example_sequence) |>
  dplyr::select(UID, Gene, Guide)


# Count how many times each guide sequence occurs
sequence_multiplicity <- library_b_data |>
  dplyr::count(
    Guide,
    name = "sequence_occurrences"
  )

# Add sequence multiplicity to the complete Library B data
library_b_data <- dplyr::left_join(
  library_b_data,
  sequence_multiplicity,
  by = "Guide"
)

# Label unique and repeated sequences
library_b_data <- dplyr::mutate(
  library_b_data,
  sequence_status = ifelse(
    sequence_occurrences == 1,
    "unique_sequence",
    "repeated_multi_gene"
  )
)

# Check the classification
table(
  library_b_data$target_type,
  library_b_data$sequence_status
)


# Examine guide coverage for each targeted gene
guide_coverage_per_gene <- library_b_data |>
  dplyr::filter(target_type == "targeting") |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    total_guides = dplyr::n(),
    unique_sequence_guides = sum(
      sequence_status == "unique_sequence"
    ),
    repeated_sequence_guides = sum(
      sequence_status == "repeated_multi_gene"
    ),
    .groups = "drop"
  )

# Current number of guides per gene
table(guide_coverage_per_gene$total_guides)

# nr that remain after excluding repeated sequences
table(guide_coverage_per_gene$unique_sequence_guides)


# Recreate subsets from the latest complete dataset
library_b_targeting <- library_b_data |>
  dplyr::filter(target_type == "targeting")


library_b_controls <- library_b_data |>
  dplyr::filter(target_type == "non_targeting_control")

# Final integrity checks
nrow(library_b_data)
nrow(library_b_targeting)
nrow(library_b_controls)

table(
  library_b_data$target_type,
  library_b_data$sequence_status
)

# Save the complete prepared Library B dataset
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  library_b_data,
  "data/processed/library_b_prepared.rds"
)

message(
  "Prepared Library B dataset saved to ",
  "data/processed/library_b_prepared.rds"
)

  
