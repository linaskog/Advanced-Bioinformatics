# Convert mapped sgRNAs to Ensembl gene IDs
# Reference genome: hg19 / GRCh37
library(EnsDb.Hsapiens.v75)
# Load our existing genome annotation results
mapped_guides<- read.delim(
  "data/processed/guide_hg19_annotated.tsv"
)
head(mapped_guides,3)

sum(complete.cases(
  mapped_guides[, c("chromosome", "start_1based", "end_1based")]
))


# Keep guides with complete genomic coordinates
mapped_guides_complete <- mapped_guides[
  complete.cases(mapped_guides[, c(
    "chromosome", "start_1based", "end_1based"
  )]),
]

# Create genomic ranges for all mapped guides
all_guide_ranges <- GenomicRanges::GRanges(
  seqnames = sub("^chr", "", mapped_guides_complete$chromosome),
  ranges = IRanges::IRanges(
    start = mapped_guides_complete$start_1based,
    end = mapped_guides_complete$end_1based
  ),
  UID = mapped_guides_complete$UID
)

length(all_guide_ranges)
# Load Ensembl gene coordinates and identifiers (GRCh37)
ensembl_genes <- ensembldb::genes(EnsDb.Hsapiens.v75)
overlap_hits<- GenomicRanges::findOverlaps(
  all_guide_ranges,
  ensembl_genes,
  ignore.strand=TRUE
  
)
length(overlap_hits)


guide_gene_ids <- data.frame(
  UID = all_guide_ranges$UID[
    queryHits(overlap_hits)
  ],
  Ensembl_ID = ensembl_genes$gene_id[
    subjectHits(overlap_hits)
  ],
  Gene_symbol = ensembl_genes$gene_name[
    subjectHits(overlap_hits)
  ]
)

head(guide_gene_ids, 10)
length(unique(guide_gene_ids$UID))
length(setdiff(all_guide_ranges$UID, guide_gene_ids$UID))

guide_id_annotations <- merge(
  data.frame(UID = all_guide_ranges$UID),
  guide_gene_ids,
  by = "UID",
  all.x = TRUE
)

head(guide_id_annotations, 10)

sum(is.na(guide_id_annotations$Ensembl_ID))

saveRDS(
  guide_id_annotations,
  "data/processed/guide_ensembl_annotations.rds"
)
