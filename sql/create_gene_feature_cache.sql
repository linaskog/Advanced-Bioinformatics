
-- PostgreSQL cache for precomputed CRISPR sgRNA features
-- Advanced Bioinformatics project
--
-- Stores feature tables associated with Ensembl gene IDs.
-- Cached features can be retrieved by Shiny without
-- repeating expensive Bowtie and genomic calculations.

CREATE TABLE IF NOT EXISTS gene_feature_cache (
    ensembl_id TEXT NOT NULL,
    gene_symbol TEXT,
    features JSONB NOT NULL,
    cache_version TEXT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,

    PRIMARY KEY (ensembl_id, cache_version)
);
