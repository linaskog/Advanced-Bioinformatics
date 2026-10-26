# Retrieve previously calculated sgRNA features from PostgreSQL

read_cached_variables <- function(con, geneID) {
  
  cache <- DBI::dbGetQuery(
    con,
    "SELECT gene_symbol, features::text AS features
     FROM gene_feature_cache
     WHERE ensembl_id = $1 AND cache_version = $2",
    params = list(geneID, "prototype_v1_unvalidated")
  )
  
  # Return NULL if the gene has not been cached
  if (nrow(cache) == 0) {
    return(NULL)
  }
  
  # Convert stored JSON back to an R data frame
  guides <- jsonlite::fromJSON(cache$features[[1]])
  
  # Restore the column order required by the models
  expected_columns <- c(
    "UID", "Guide",
    paste0(
      "pos",
      rep(1:20, each = 4),
      "_",
      rep(c("A", "C", "G", "T"), 20)
    ),
    "GC", "strand", "chromosome",
    "region_type", "gene_overlap_count",
    "H3K27ac_peak_overlap"
  )
  
  if (!setequal(names(guides), expected_columns)) {
    stop("Cached features do not match the expected model inputs.")
  }
  
  guides <- guides[, expected_columns, drop = FALSE]
  
  # Restore gene metadata
  attr(guides, "gene_symbol") <- cache$gene_symbol[[1]]
  attr(guides, "ensembl_id") <- geneID
  
  return(guides)
}


# Retrieve cached features or calculate and save them

get_cached_variables <- function(con, geneID) {
  
  # Step 1: Check PostgreSQL
  guides <- read_cached_variables(con, geneID)
  
  if (!is.null(guides)) {
    message("Cache HIT: ", geneID)
    return(guides)
  }
  
  # Step 2: Calculate features if not cached
  message("Cache MISS: calculating features for ", geneID)
  
  guides <- get_variables(geneID)
  
  if (!is.data.frame(guides) || nrow(guides) == 0) {
    stop("No valid guide features available to cache.")
  }
  
  # Step 3: Convert results to JSON
  features_json <- jsonlite::toJSON(
    guides,
    dataframe = "rows",
    na = "null",
    auto_unbox = TRUE,
    digits = NA
  )
  
  # Step 4: Save results in PostgreSQL
  rows_inserted<- DBI::dbExecute(
    con,
    "INSERT INTO gene_feature_cache
     (ensembl_id, gene_symbol, features, cache_version)
     VALUES ($1, $2, $3::jsonb, $4)
     ON CONFLICT (ensembl_id, cache_version) DO NOTHING",
    params = list(
      geneID,
      attr(guides, "gene_symbol"),
      as.character(features_json),
      "prototype_v1_unvalidated"
    )
  )
  
  
  if (rows_inserted == 1L) {
    message("Features saved to PostgreSQL: ", geneID)
  } else {
    message("Cache entry already exists; insert skipped: ", geneID)
  }
  
  
  return(guides)
}