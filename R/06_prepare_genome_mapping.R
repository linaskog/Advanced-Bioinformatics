#read library_b_prepared.rds
#select targeting guides
#verify 20-nt A/C/G/T sequences
#create FASTA
#write library_b_targeting_guides.fa

library_b<- readRDS("data/processed/library_b_prepared.rds")
head(library_b)
dim(library_b)
names(library_b)
head(library_b[, c("UID", "Guide", "Gene", "target_type")])
table(library_b$target_type)
targeting_guides<- library_b[library_b$target_type=="targeting",]
dim(targeting_guides)
head(targeting_guides)

table(nchar(targeting_guides$Guide))
table(grep("^[ACGT]+$", targeting_guides$Gene))
table(grepl("^[ACGT]+$", targeting_guides$Guide))

# fasta file (targeting_guides)
fasta_records<- paste0(">", targeting_guides$UID,
                      "\n",
                      targeting_guides$Guide)

head(fasta_records)
length(fasta_records)

#save the fasta_record 
writeLines(
  fasta_records,
  "data/processed/library_b_targeting_guides.fa"
)


