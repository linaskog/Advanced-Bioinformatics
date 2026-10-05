library(rtracklayer)
library(GenomicRanges)
library(IRanges)

gtf<- rtracklayer::import(
  "data/reference/hg19.refGene.gtf.gz"
)
class(gtf)
mcols(gtf)
colnames(mcols(gtf))
transcripts<- gtf[gtf$type== "transcript"]
exons <- gtf[gtf$type == "exon"]
length(transcripts)
length(exons)

guide_coords<- read.delim(
  "data/processed/guide_hg19_unique_coordinates.tsv",
  stringsAsFactors = FALSE
)
dim(guide_coords)
colnames(guide_coords)

#convert data frame into a GRanges object,
#so R can treat each sgRNA as a genomic interval.
guide_gr <- GRanges(
  seqnames = guide_coords$chromosome,
  ranges = IRanges(
    start = guide_coords$start_1based,
    end = guide_coords$end_1based 
  ),
  strand = guide_coords$strand,
  UID= guide_coords$UID
)

length(guide_gr)

# Region classification:
# exon = overlaps at least one annotated exon
# intron = overlaps a transcript but no annotated exon
# intergenic = overlaps no annotated transcript
guide_transcript_hits <- findOverlaps(
  guide_gr,
  transcripts,
  ignore.strand=TRUE
)
length(guide_transcript_hits)
length(unique(queryHits(guide_transcript_hits)))

#guides overlap exons
guide_exon_hits <- findOverlaps(
  guide_gr,
  exons,
  ignore.strand=TRUE
)

length(guide_exon_hits)
length(unique(queryHits(guide_exon_hits)))


#calculate introns
transcript_guides <- unique(queryHits(guide_transcript_hits))
exon_guides <- unique(queryHits(guide_exon_hits))
intron_only_guides <- setdiff(transcript_guides, exon_guides)
length(intron_only_guides)

#intergenic guides
intergenic_guides<- setdiff(
  seq_along(guide_gr),
  transcript_guides
)
length(intergenic_guides)

# add region_type clumn to the  guide_coords
guide_coords$region_type <- "intergenic"

guide_coords$region_type[intron_only_guides] <- "intron"


guide_coords$region_type[exon_guides] <- "exon"
table(guide_coords$region_type)

#Which gene does each guide overlap
guide_gene_pairs <- data.frame(
  UID = guide_gr$UID[queryHits(guide_transcript_hits)],
  gene_name= transcripts$gene_name[subjectHits(guide_transcript_hits)],
  stringsAsFactors = FALSE
)

head(guide_gene_pairs)
guide_gene_pairs_unique <- unique(guide_gene_pairs)
head(guide_gene_pairs_unique)
nrow(guide_gene_pairs_unique)

genes_per_guide<- table(guide_gene_pairs_unique$UID)
table(genes_per_guide)
sum(genes_per_guide == 1)
sum(genes_per_guide > 1)

guide_gene_collapsed <- aggregate(
  gene_name ~ UID,
  data = guide_gene_pairs_unique,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)
 
head(guide_gene_collapsed)

gene_overlap_count <- as.integer(genes_per_guide[
  match(guide_coords$UID, names(genes_per_guide))
])

guide_coords$gene_overlap_count <- gene_overlap_count

guide_coords$gene_overlap_count[
  is.na(guide_coords$gene_overlap_count)
] <- 0

table(guide_coords$gene_overlap_count)

# attaching gene annotation to the table
guide_coords$gene_name <- guide_gene_collapsed$gene_name[
  match(guide_coords$UID, guide_gene_collapsed$UID)
]

head(guide_coords)
sum(is.na(guide_coords$gene_name))

stopifnot(nrow(guide_coords) == 52357)

stopifnot(
  sum(guide_coords$region_type == "exon") == 52281,
  sum(guide_coords$region_type == "intron") == 26,
  sum(guide_coords$region_type == "intergenic") == 50
)

stopifnot(sum(guide_coords$gene_overlap_count == 0) == 50)

write.table(
  guide_coords,
  "data/processed/guide_hg19_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
