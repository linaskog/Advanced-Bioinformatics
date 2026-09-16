# Advanced Bioinformatics Project

## Project overview

This project analyses data from a genome-wide CRISPR-Cas9 knockout screen performed in two B-cell cancer cell lines:

* **BC-3:** a KSHV-positive primary effusion lymphoma cell line.
* **BJAB:** a KSHV-negative Burkitt lymphoma cell line.

The screen used the **GeCKO v2 sgRNA library**. The library is divided into Library A and Library B. Together, they generally contain six sgRNAs per protein-coding gene, with approximately three guides from each library.

The long-term objective is to develop a model that predicts sgRNA performance or quality from sgRNA sequence features. The precise response variable used to represent sgRNA quality has not yet been finalized.

## Biological principle

CRISPR-Cas9 is used to knock out genes in a population of cells. The abundance of each sgRNA is measured at the beginning of the experiment, Day 0, and after selection, Day 14.

If an sgRNA becomes less abundant at Day 14, disrupting its target may have reduced cell growth or survival. However, the observed change can be influenced by:

* sgRNA efficiency;
* biological importance of the target gene;
* cell-line-specific effects;
* sequencing depth;
* experimental noise.

Therefore, raw fold change cannot automatically be interpreted as sgRNA quality. The response variable and validation strategy must be selected carefully.

## Data source

The data originate from the supplementary material of Manzano et al. (2018):

* `41467_2018_5506_MOESM4_ESM.xlsx`
* `41467_2018_5506_MOESM6_ESM.xlsx`

Associated publication:

[Manzano et al., Nature Communications, 2018](https://doi.org/10.1038/s41467-018-05506-9)

The raw Excel files are stored locally in:

```text
data/raw/
```

The raw data are excluded from Git using `.gitignore`. They must not be modified directly.

## Imported datasets

The script `R/01_import_data.R` imports four datasets.

| R object        | Excel sheet         | Contents                                          |   Dimensions |
| --------------- | ------------------- | ------------------------------------------------- | -----------: |
| `gecko_counts`  | `C. GECKO_BC3_BJAB` | sgRNA counts for BC-3 and BJAB                    | 123,411 × 14 |
| `gecko_library` | `F. GeCKOlibrary2`  | sgRNA identifiers, sequences and gene annotations |  123,411 × 3 |
| `bc3_mageck`    | `GECKO_BC-3`        | Gene-level MAGeCK results for BC-3                |  21,915 × 14 |
| `bjab_mageck`   | `GECKO_BJAB`        | Gene-level MAGeCK results for BJAB                |  21,915 × 14 |

The count table contains:

* BC-3 Day 0: three replicates;
* BC-3 Day 14: three replicates;
* BJAB Day 0: three replicates;
* BJAB Day 14: three replicates.

## Project scripts

```text
R/
├── 00_exploration.R
├── 01_import_data.R
└── 02_quality_control.R
```

### `00_exploration.R`

Contains the initial exploratory work used to understand the Excel files and identify the relevant worksheets.

### `01_import_data.R`

Imports the required worksheets and verifies that the expected columns are present.

### `02_quality_control.R`

Performs quality-control checks on the imported data, including:

1. dataset dimensions;
2. column data types;
3. missing values;
4. duplicated identifiers;
5. agreement of identifiers between tables;
6. agreement of gene annotations;
7. negative and zero counts;
8. sequencing depth;
9. agreement between biological replicates.

## Quality-control results

### Dimensions and data types

All four datasets were imported with the expected dimensions.

* UID, gene and guide-sequence columns were imported as character data.
* Experimental counts and MAGeCK statistics were imported as numeric data.
* BC-3 and BJAB MAGeCK tables had identical column structures.

### Missing values

No missing values were detected in any of the four imported datasets.

### Duplicate identifiers

No duplicated identifiers were found in:

* `gecko_counts$UID`;
* `gecko_library$UID`;
* `bc3_mageck$id`;
* `bjab_mageck$id`.

### UID coverage

Both sgRNA-level tables contained the same set of 123,411 UIDs:

* count UIDs missing from the sequence library: 0;
* library UIDs missing from the count table: 0.

Although the UIDs were present in both tables, their associated gene annotations were not always consistent.

## UID-to-gene mapping problem

Matching the count and sequence tables by UID initially produced:

| Mapping result            | Number of sgRNAs | Percentage |
| ------------------------- | ---------------: | ---------: |
| Matching gene annotation  |           81,175 |     65.78% |
| Different gene annotation |           42,236 |     34.22% |

The mismatch was systematic rather than random.

### Library-specific results

| GeCKO library | Total guides | Matching guides | Mismatching guides | Match rate |
| ------------- | -----------: | --------------: | -----------------: | ---------: |
| Library A     |       65,383 |          23,150 |             42,233 |     35.41% |
| Library B     |       58,028 |          58,025 |                  3 |     99.99% |

Almost all UID-to-gene inconsistencies occurred in Library A.

The three apparent Library B mismatches involved the gene `DEC1`. Excel represented the same value differently:

* `1-Dec` in the count table;
* `42339` in the sequence-library table.

After correcting this formatting problem in derived variables, Library B had zero remaining UID-to-gene mismatches. The original raw data were not modified.

Library B can therefore be connected to the guide sequences by UID. Library A must not be connected to guide sequences until its mapping problem has been resolved or an independently verified mapping source has been obtained.

## miRNA targets

The dataset contains 1,864 miRNA targets. These were represented by 7,288 mismatched sgRNA rows.

If every mismatched guide were excluded:

* 42,236 sgRNAs would be removed;
* 81,175 sgRNAs would remain;
* 1,867 targets would disappear completely;
* the lost targets would include all 1,864 miRNA targets, plus `DEC1`, `PLN` and `SPHAR`.

The original publication did not report the miRNA-screen results because miRNA knockout analysis has additional biological complications.

The miRNA mismatches are part of the Library A mapping problem. They are not an independent source of mismatch.

## Count quality

No negative counts were detected in any experimental sample.

Zero counts were uncommon in most samples. BJAB Day 14 contained more zero-count guides:

| Experimental group | Percentage of zero counts |
| ------------------ | ------------------------: |
| BC-3 Day 0         |                0.01–0.02% |
| BC-3 Day 14        |                0.22–0.29% |
| BJAB Day 0         |                0.07–0.14% |
| BJAB Day 14        |                3.34–6.85% |

`BJAB_Day14_Rep1` had the largest percentage of zero-count guides: 6.85%.

Zero counts are not automatically errors. Some guides may disappear during selection because their target genes affect cell survival. However, differences between replicates must be considered during subsequent analysis.

## Sequencing depth

Total sequencing depth ranged from approximately 9.83 million to 29.24 million reads per sample.

These differences show that normalization will be required before comparing counts or calculating sgRNA-level changes. No sample was excluded based only on sequencing depth.

## Replicate agreement

Pearson correlations were calculated using log2-transformed counts:

$$
\log_2(\text{count}+1)
$$

The within-group replicate correlations were:

| Experimental group | Correlation range | Interpretation             |
| ------------------ | ----------------: | -------------------------- |
| BC-3 Day 0         |         0.94–0.96 | Very strong agreement      |
| BC-3 Day 14        |         0.80–0.83 | Moderate-to-good agreement |
| BJAB Day 0         |         0.88–0.89 | Strong agreement           |
| BJAB Day 14        |         0.54–0.58 | Weak agreement             |

BJAB Day 14 showed the weakest replicate agreement and the largest proportion of zero counts. These samples will not be removed automatically, but their variability must be considered when defining the model response and evaluating model performance.

## Current conclusions

The following conclusions have been reached:

1. The required datasets can be imported reproducibly.
2. No missing values, duplicate identifiers or negative counts were detected.
3. Sequencing depth differs between samples, so normalization is required.
4. BC-3 and BJAB Day 0 replicates show strong agreement.
5. BJAB Day 14 replicates show relatively weak agreement.
6. Library B has a reliable UID-to-sequence mapping after correction of the `DEC1` Excel-formatting issue.
7. Library A contains a systematic UID-to-gene mapping inconsistency and cannot yet be safely used for guide-level sequence modeling.
8. The dataset includes 1,864 miRNA targets, which require separate consideration.
9. No raw data have been modified or deleted.


## Reproducing the current analysis

Open the R project and ensure that the raw Excel files are available in `data/raw/`.

Install the required import package if necessary:

```r
install.packages("readxl")
```

Run the quality-control analysis from the project root:

```r
source("R/02_quality_control.R")
```

The quality-control script sources `R/01_import_data.R`, so the required data are imported before the checks are performed.

## Project status

The import and initial quality-control stages are complete.

No final modeling dataset or prediction model has been created yet. The next decision is how to handle Library A and whether the first modeling analysis should use only the verified Library B guides.

