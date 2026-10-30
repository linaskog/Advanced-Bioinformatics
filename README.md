# Advanced Bioinformatics — CRISPR Screen Analysis and sgRNA Scoring

A course project developed for **Advanced Bioinformatics (7.5 ECTS), Umeå University**.

## Overview

We analysed a genome-wide CRISPR-Cas9 knockout screen in two B-cell lymphoma cell lines: **BC-3** (KSHV-positive primary effusion lymphoma) and **BJAB** (KSHV-negative Burkitt lymphoma). The project combines screening-data quality control, genomic annotation, database development and machine learning in an interactive R Shiny application.

## Project aim

The aim of this project is to develop a bioinformatics tool for identifying, evaluating, and ranking candidate sgRNAs for CRISPR-Cas9 gene editing, inspired by tools such as CHOPCHOP.

Given an Ensembl gene ID, the application identifies potential sgRNA sequences and calculates sequence- and genome-related features, including GC content, genomic location, and potential off-target matches. Linear regression and XGBoost models are used to predict scores and rank candidate guides.

The tool is implemented as an interactive R Shiny application, allowing users to explore candidate sgRNAs, assess sequence quality, and compare model predictions.

## Data preparation and quality control

The screening data come from the supplementary material of [Manzano et al. (2018)](https://doi.org/10.1038/s41467-018-05506-9). The GeCKO v2 dataset contains **123,411 sgRNAs**, with counts from three replicates at Day 0 and Day 14 for each cell line.

We checked sample quality, sequencing depth, missing values, gene annotations and agreement between replicates. A substantial UID-to-gene annotation mismatch was found in **Library A**, whereas Library B was consistent after correcting an Excel formatting issue affecting `DEC1`. We therefore used **Library B** for the downstream guide-level workflow.

The prepared Library B dataset contains **58,028 guides**: **57,028 targeting sgRNAs** and **1,000 non-targeting controls**. We retained controls and labelled repeated guide sequences rather than removing them automatically. Sample-level QC also identified weaker replicate agreement in **BJAB Day 14** (approximately 0.54–0.58), which is important when interpreting screening results.

- Preparation: `R/03_prepare_analysis_data.R`
- Sample QC: `R/04_sample_level_qc.R`
- Prepared data: `data/processed/library_b_prepared.rds`
- QC figures and summary: `results/qc/`

## Genome mapping and annotation

To add genomic context, we mapped the **57,028 targeting guides** to the **hg19/GRCh37** reference genome with **Bowtie 1.3.1**, using **HPC2N/SLURM**. Alignments required an exact sequence match, and all matching genomic locations were reported. Of these guides, **52,357 mapped uniquely** and were carried forward for location-based annotation.

Using the hg19 RefGene annotation and genomic overlap operations in R, we classified unique guide locations as **exonic, intronic or intergenic** and recorded overlaps with annotated genes. We also integrated **BC-3 H3K27ac ChIP-seq peaks** as a binary genomic feature; **4,118 uniquely mapped guides** overlapped a peak. Additional Ensembl-based annotation was explored to connect mapped coordinates with gene IDs and symbols.

Relevant scripts are `R/06_prepare_genome_mapping.R`, `R/07_annotate_genome_mapping.R` and `R/08_integrate_bc3_h3k27ac.R`. The mapping jobs are in `slurm/`, and processed annotations are in `data/processed/`.

## Databases

### SQLite: screening data and QC

We built a relational SQLite database to keep guide information, sample metadata, raw counts and sample-level QC in one structured resource. Its four connected tables are `guides`, `samples`, `counts` and `sample_qc`. Together they contain **58,028 guides**, **12 samples** and **696,336 guide-by-sample count records**. Primary and foreign keys, data constraints and validation checks help maintain consistency.

- Build script: `R/05_build_sqlite_database.R`
- Database: `data/processed/crispr_screen.sqlite`

### PostgreSQL: Shiny feature cache

The Shiny application uses PostgreSQL to reuse genomic features calculated for a selected Ensembl gene ID. On a **cache hit**, previously calculated features are loaded; on a **cache miss**, the application calculates and stores them for later use. The implementation uses parameterized SQL queries and stores feature tables in a JSONB column, with a cache version to distinguish calculation versions.

- Cache functions: `R/app_postgres_cache.R`
- Table definition: `sql/create_gene_feature_cache.sql`
- Configuration example: `.Renviron.example`

The PostgreSQL cache was tested locally, including successful retrieval of previously computed GATA3 features through Shiny. Cache misses trigger the original feature-calculation pipeline; this route is still under validation, particularly for guides with extensive genomic multi-mapping. Cached features use the version label `prototype_v1_unvalidated`.

## Statistical and Machine Learning Models and Evaluation

We developed two predictive models, **linear regression and XGBoost**, to estimate sgRNA performance using CRISPR screening data from BC-3 cells. Both models were trained to predict gene-centered log2 fold change using sgRNA sequence information, GC content, and genomic features such as chromosome, strand, region type, gene overlaps, and H3K27ac peak overlaps.

The models were evaluated using a separate test dataset. The best-performing linear regression model achieved a test R² of approximately **0.046**, while XGBoost achieved **0.050**. These results show that both models have limited predictive performance, although XGBoost performed slightly better.

For more information about model development, training, hyperparameter optimization, evaluation, and SHAP analysis, see the [detailed modelling documentation](models.md).

The model scripts are available in `R/`, the saved models in `models/`, and the evaluation figures in `results/models/`. The repository also contains a small `R/basicstats/` package for model evaluation metrics. 


## Shiny application

The primary R Shiny interface (`R/test_app.R`) provides three tabs:

1. **Input Gene ID** — Enter an Ensembl gene ID to generate candidate 20-nucleotide sgRNAs.
2. **Quality Control** — Review guide counts, mean GC content, an adjustable GC-content distribution and an interactive table. The table flags sequences outside the selected GC range, poly-T (`TTTT`) motifs and other homopolymer runs of four or more bases.
3. **Score** — Compare linear-model and XGBoost predictions, adjust their relative weights, rank candidate guides and inspect an interactive plot and results table.

Sequence-quality flags are screening aids, **not direct measurements of editing efficiency or off-target specificity**. The application retrieves feature inputs through the PostgreSQL cache. The earlier `R/new_app.R` script remains in the repository, but `R/test_app.R` is the current integrated application.

## Getting started

Open `Advanced-Bioinformatics.Rproj` in RStudio and use the repository root as the working directory.

**Requirements:** R and relevant CRAN/Bioconductor packages (including `shiny`, `bslib`, `DBI`, `RSQLite`, `RPostgres`, `jsonlite`, `xgboost`, `DT`, `plotly`, `GenomicRanges`, `rtracklayer`, `EnsDb.Hsapiens.v75`, `BSgenome.Hsapiens.UCSC.hg19` and `Rbowtie`). The full genome-mapping workflow also requires Bowtie and an hg19 reference/index.

**To reproduce data preparation**, place the original supplementary Excel files in `data/raw/` and run the relevant scripts from the project root, for example:

```r
source("R/02_quality_control.R")
source("R/03_prepare_analysis_data.R")
source("R/04_sample_level_qc.R")
source("R/05_build_sqlite_database.R")
```

**To run the Shiny app locally:**

1. Set up a PostgreSQL database and create the cache table using `sql/create_gene_feature_cache.sql`.
2. Use `.Renviron.example` to configure `CRISPR_DB_HOST`, `CRISPR_DB_PORT`, `CRISPR_DB_NAME`, `CRISPR_DB_USER` and `CRISPR_DB_PASSWORD` in a local `.Renviron` file. Do not commit credentials or overwrite an existing local configuration.
3. Make sure the hg19 Bowtie index (prefix `hg19_index/hg19`), reference annotations, and saved models are available. The active Shiny scripts now use paths relative to the **repository root**; some separate model-development scripts still contain computer-specific paths.
4. Launch Shiny from the RStudio Console while the working directory is the repository root:

```r
app <- source("R/test_app.R", chdir = FALSE)$value
shiny::runApp(app)
```

Launching `test_app.R` directly with RStudio’s **Run App** button may change the working directory to `R/`, which breaks paths such as `source("R/app_get_variables.R")`. This setup has been tested locally; it is **not yet a fully automated installation procedure**.

## Repository guide

| Location | Contents |
| --- | --- |
| `R/` | Import, QC, mapping, modeling, database and Shiny scripts |
| `R/basicstats/` | Custom model-evaluation package |
| `data/processed/` | Prepared guides, annotations and SQLite database |
| `data/reference/` | Genome annotation files |
| `models/` | Saved linear and XGBoost models |
| `results/qc/`, `results/models/` | QC and model figures |
| `slurm/` | HPC genome-indexing and guide-mapping jobs |
| `sql/` | PostgreSQL cache-table definition |

## Current limitations

Genome-based feature generation on a cache miss still needs validation. In particular, highly multi-mapping guides can produce very large Bowtie outputs, and guide-to-gene overlap counting and genomic coordinate conventions are being reviewed in a separate experimental script. These experimental fixes are **not** part of the integrated application. Existing cached features should not be interpreted as fully validated. PostgreSQL is local to each installation, and required genomic assets must be provided separately. The modeling contributor should add the final evaluation metrics before submission.

## Reference

Manzano, M., Patil, A., Waldrop, A., Dave, S. S., Behdad, A., & Gottwein, E. (2018). Gene essentiality landscape and druggable oncogenic dependencies in herpesviral primary effusion lymphoma. Nature Communications, 9, Article 3263. https://doi.org/10.1038/s41467-018-05506-9

Manzano, M., Patil, A., Waldrop, A., Dave, S. S., Behdad, A., & Gottwein, E. (2018). Supplemental Data Legends [Supplementary file accompanying the article above]. Nature Communications. https://doi.org/10.1038/s41467-018-05506-9
