
library(rtracklayer)
library(GenomicRanges)

h3k27ac_peaks <- rtracklayer::import(
  "data/raw/GSM4040945_BC3.merge.peak.bed.gz",
  format = "BED"
)
head(h3k27ac_peaks)
class(h3k27ac_peaks)
length(h3k27ac_peaks)

guide_coords <- read.delim(
  "data/processed/guide_hg19_annotated.tsv",
  stringsAsFactors = FALSE
)
head(guide_coords)

guide_gr <- GRanges(
  seqnames = guide_coords$chromosome,
  ranges = IRanges(
    start = guide_coords$start_1based,
    end = guide_coords$end_1based
  ),
  strand = guide_coords$strand,
  UID = guide_coords$UID
)
length(guide_gr)

# find the overlaps region
guide_h3k27ac_hits <- findOverlaps(
  guide_gr,
  h3k27ac_peaks,
  ignore.strand = TRUE
)
length(guide_h3k27ac_hits)
length(unique(queryHits(guide_h3k27ac_hits)))


guide_coords$H3K27ac_peak_overlap <- FALSE

guide_coords$H3K27ac_peak_overlap[
  unique(queryHits(guide_h3k27ac_hits))
] <- TRUE

table(guide_coords$H3K27ac_peak_overlap)

stopifnot(nrow(guide_coords) == 52357)

stopifnot(length(h3k27ac_peaks) == 25546)

stopifnot(
  sum(guide_coords$H3K27ac_peak_overlap) == 4118
)

write.table(
  guide_coords,
  "data/processed/guide_hg19_bc3_h3k27ac.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
