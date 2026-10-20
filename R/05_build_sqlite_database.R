# Build the CRISPR-screen SQLite database

library(DBI)
library(RSQLite)

prepared_data_file <- "data/processed/library_b_prepared.rds"
sample_qc_file <- "results/qc/library_b_sample_qc_summary.csv"

stopifnot(
  file.exists(prepared_data_file),
  file.exists(sample_qc_file)
)

library_b_data <- readRDS(
  prepared_data_file
)

sample_qc <- read.csv(
  sample_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(library_b_data)
dim(sample_qc)
names(sample_qc)

# Identify the 12 experimental count columns
count_columns <- grep(
  "^(BC-3|BJAB)_Day(0|14)_Rep[123]$",
  names(library_b_data),
  value = TRUE
)

# Confirm that all 12 samples were found and match the QC table
stopifnot(
  length(count_columns) == 12,
  setequal(count_columns, sample_qc$sample)
)



#samples table
samples_table <- data.frame(
  sample_id= seq_along(count_columns),
  sample= count_columns,
  cell_line= sub(
    "_Day.*$",
    "",
    count_columns
  ),
  day= sub(
    "^.*_(Day0|Day14)_.*$",
    "\\1",
    count_columns
  ),
  replicate= sub(
    ".*_(Rep[123])$",
    "\\1",
    count_columns
  )
)

stopifnot(
  all(
    table(samples_table$cell_line) ==
      c("BC-3" = 6L, "BJAB" = 6L)
  ),
  all(
    table(samples_table$day) ==
      c("Day0" = 6L, "Day14" = 6L)
  )
)

samples_table


#guides table
guides_table<- data.frame(
  uid= library_b_data$UID,
  guide_sequence= library_b_data$Guide,
  gene=library_b_data$Gene,
  target_type=library_b_data$target_type,
  sequence_occurrences= library_b_data$sequence_occurrences,
  sequence_status=library_b_data$sequence_status
)

stopifnot(
  nrow(guides_table)==58028,
  !anyDuplicated(guides_table$uid),
  !anyNA(guides_table)
)

dim(guides_table)
head(guides_table)

# sample_qc table
#this table connect the qc_recored to the sample id
sample_qc_table <- data.frame(
  sample_id= samples_table$sample_id[
    match(sample_qc$sample, samples_table$sample)
  ],
  total_reads= sample_qc$total_reads,
  zero_count_percentage= sample_qc$zero_count_percentage,
  deseq2_size_factor= sample_qc$deseq2_size_factor,
  PC1= sample_qc$PC1,
  PC2= sample_qc$PC2
)

stopifnot(
  nrow(sample_qc_table)==12,
  !anyDuplicated(sample_qc_table$sample_id),
  !anyNA(sample_qc_table)
)

sample_qc_table


# counts table
count_matrix<- as.matrix(
  library_b_data[, count_columns]
)

counts_table <- data.frame(
  uid= rep(
    library_b_data$UID,
    times=length(count_columns)
  ),
  sample_id= rep(
    samples_table$sample_id,
    each=nrow(library_b_data)
    ),
    raw_count= as.integer(count_matrix)
)

stopifnot(
  nrow(counts_table)== 58028*12,
  !anyNA(counts_table)
)
count_totals_check <- aggregate(
  raw_count ~ sample_id,
  data= counts_table,
  FUN = sum
)

count_totals_check
dim(counts_table)

# Create the SQLite database
stopifnot(
  all(
    count_totals_check$raw_count==
      sample_qc_table$total_reads
  )
)

database_file <- "data/processed/crispr_screen.sqlite"
if(file.exists(database_file)) {
  file.remove(database_file)
}

connection<- dbConnect(       #connection between R and database
  SQLite(),
  database_file
)


dbExecute(
  connection,
  "PRAGMA foreign_keys =ON;"
)

dbGetQuery(
  connection,
  "PRAGMA foreign_keys;"
)

# Create and populate the samples table
dbExecute(
  connection,
  "
  CREATE TABLE samples(
    sample_id INTEGER PRIMARY KEY,
    sample TEXT NOT NULL UNIQUE,
    cell_line TEXT NOT NULL
      CHECK(cell_line IN ('BC-3', 'BJAB')),
    day TEXT NOT NULL
      CHECK(day IN ('Day0', 'Day14')),
    replicate TEXT NOT NULL
      CHECK (replicate IN ('Rep1', 'Rep2', 'Rep3'))
  );
  "
)

dbAppendTable(
  connection,
  "samples",
  samples_table
)

dbGetQuery(
  connection,
  "SELECT * FROM samples"
)


# creat table guide
dbExecute(
  connection,
  "
  CREATE TABLE guides (
    uid TEXT PRIMARY KEY,
    guide_sequence TEXT NOT NULL
      CHECK (length(guide_sequence) = 20),
    gene TEXT NOT NULL,
    target_type TEXT NOT NULL
      CHECK (
        target_type IN (
          'targeting',
          'non_targeting_control'
        )
      ),
    sequence_occurrences INTEGER NOT NULL
      CHECK (sequence_occurrences >= 1),
    sequence_status TEXT NOT NULL
      CHECK (
        sequence_status IN (
          'unique_sequence',
          'repeated_multi_gene'
        )
      )
  );
  "
)
dbAppendTable(
  connection,
  "guides",
  guides_table
)

dbGetQuery(
  connection,
  "
  SELECT
    target_type,
    sequence_status,
    COUNT(*) AS number_of_guides
  FROM guides
  GROUP BY target_type, sequence_status;
  "
)


#sample-QC table
dbExecute(
  connection,
  "
  CREATE TABLE sample_qc(
  sample_id INTEGER PRIMARY KEY,
  total_reads INTEGER NOT NULL
    CHECK(total_reads >= 0),
  zero_count_percentage REAL NOT NULL
    CHECK(
      zero_count_percentage >= 0
      AND zero_count_percentage <= 100
    ),
  deseq2_size_factor REAL NOT NULL
    CHECK(deseq2_size_factor > 0),
  PC1 REAL NOT NULL,
  PC2 REAL NOT NULL,
  FOREIGN KEY (sample_id)
  REFERENCES samples(sample_id)
  ON DELETE CASCADE
  
  );
  "
)

dbAppendTable(
  connection,
  "sample_qc",
  sample_qc_table
)
dbGetQuery(
  connection,
  "
  SELECT
    samples.sample,
    samples.cell_line,
    samples.day,
    sample_qc.zero_count_percentage
  FROM samples
  INNER JOIN sample_qc
    ON samples.sample_id= sample_qc.sample_id
  ORDER BY samples.sample_id;  
  "
)

# counts table (the largest table  696,336)
dbExecute(
  connection,
  "
  CREATE TABLE counts(
  uid TEXT NOT NULL,
  sample_id INTEGER NOT NULL,
  raw_count INTEGER NOT NULL
    CHECK(raw_count>= 0),
    PRIMARY KEY (uid, sample_id),
    FOREIGN KEY (uid)
      REFERENCES guides(uid)
      ON DELETE CASCADE,
    FOREIGN KEY (sample_id)
      REFERENCES samples(sample_id)
      ON DELETE CASCADE
  ) WITHOUT ROWID;
  "
)

dbAppendTable(
  connection,
  "counts",
  counts_table
)

dbExecute(
  connection,
  "
  CREATE INDEX counts_sample_id_index
  ON counts(sample_id);
  "
)

dbGetQuery(
  connection,
  "
  SELECT
    samples.sample,
    SUM(counts.raw_count) AS total_reads
  FROM counts
  INNER JOIN samples
    ON counts.sample_id = samples.sample_id
  GROUP BY samples.sample_id, samples.sample
  ORDER BY samples.sample_id
  "
)
# Validate and close the database ----------------------------------------

database_table_sizes <- dbGetQuery(
  connection,
  "
  SELECT 'guides' AS table_name, COUNT(*) AS rows
  FROM guides

  UNION ALL

  SELECT 'samples', COUNT(*)
  FROM samples

  UNION ALL

  SELECT 'sample_qc', COUNT(*)
  FROM sample_qc

  UNION ALL

  SELECT 'counts', COUNT(*)
  FROM counts;
  "
)

foreign_key_issues <- dbGetQuery(
  connection,
  "PRAGMA foreign_key_check;"
)

print(database_table_sizes)
print(foreign_key_issues)

stopifnot(
  nrow(foreign_key_issues) == 0
)

dbDisconnect(connection)

message(
  "SQLite database created successfully at ",
  database_file
)
