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

## Repository structure

```text
Advanced-Bioinformatics/
├── Advanced-Bioinformatics.Rproj
├── README.md
├── data/
│   ├── raw/
│   │   ├── 41467_2018_5506_MOESM4_ESM.xlsx
│   │   └── 41467_2018_5506_MOESM6_ESM.xlsx
│   └── processed/
│       └── library_b_prepared.rds
└── R/
    ├── 00_exploration.R
    ├── 01_import_data.R
    ├── 02_quality_control.R
    └── 03_prepare_analysis_data.R
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

### `03_prepare_analysis_data.R`

Prepares Library B for downstream analysis. The script:

1. classifies guides as Library A or Library B from their UID prefixes;
2. filters the count and sequence tables independently to Library B;
3. joins the tables safely by UID;
4. corrects the `DEC1` Excel-formatting inconsistency in derived annotation columns;
5. validates guide sequences and experimental counts;
6. classifies targeting guides and non-targeting controls;
7. identifies and labels repeated guide sequences without removing them; and
8. creates the complete cleaned Library B dataset used for subsequent QC and modelling preparation.
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

## Library selection and preparation

The UID prefix identifies the GeCKO sublibrary:

| UID prefix | Library | Number of guides |
| ---------- | ------- | ----------------: |
| `HGLibA_`  | Library A | 65,383 |
| `HGLibB_`  | Library B | 58,028 |

All 123,411 guides were assigned to one of these libraries, with no unknown UID prefixes.

Library B is currently used for downstream guide-level analysis because its count and sequence tables can be connected reliably by UID after correction of the `DEC1` formatting problem. Library A has not been permanently discarded. Its annotation structure, controls, sequences, count completeness and sample-level behaviour will be investigated separately before a final decision is made.

For Library B, the count and sequence tables were filtered independently using their own UID columns and then joined by `UID`. This avoids relying on row order or requiring the original gene annotations to agree during the join.

The joined Library B dataset contained:

- 58,028 sgRNA rows;
- no missing guide sequences;
- no missing experimental counts;
- no duplicated UIDs;
- 58,025 initially matching gene annotations; and
- three apparent `DEC1` annotation mismatches caused by Excel formatting.

Clean annotation columns, `Gene_counts_clean` and `Gene_library_clean`, were created without changing the imported raw data. After correcting `1-Dec` and `42339` to `DEC1`, all 58,028 Library B annotations agreed. The final analysis annotation is stored in the `Gene` column.

## Library B sequence validation

All 58,028 Library B guide sequences passed the sequence-level validation checks:

- every guide sequence was 20 nucleotides long;
- every sequence contained only `A`, `C`, `G` and `T`;
- missing guide sequences: 0;
- missing experimental counts: 0; and
- duplicated UIDs: 0.

These checks confirm that the cleaned Library B dataset is structurally complete for subsequent quality control and sequence-feature generation.

## Targeting guides and non-targeting controls

Library B contains:

| Guide type | Number of guides |
| ---------- | ---------------: |
| Targeting guides | 57,028 |
| Non-targeting controls | 1,000 |
| **Total** | **58,028** |

The targeting guides represent 19,049 annotated genes. The controls are identified by annotations beginning with:

```text
NonTargetingControlGuideForHuman_
```
A `target_type` column was created with two possible values:

- `targeting`
- `non_targeting_control`

The controls were retained in the prepared dataset. They may be useful for quality control, normalization, negative-reference distributions and sensitivity analyses. They will not automatically be included as ordinary gene-targeting observations during model training.

No miRNA annotations were detected within Library B using an annotation search.


## Repeated guide sequences in Library B

Guide sequences were checked independently of UID duplication. Although every UID was unique, some 20-nucleotide guide sequences occurred in multiple rows.

The repeated-sequence analysis found:

- 1,159 extra occurrences beyond the first occurrence;
- 1,903 rows involved in sequence repetition;
- 744 distinct repeated sequences;
- no repeated sequences among the 1,000 non-targeting controls; and
- repeated sequences assigned to between 2 and 21 gene annotations.

All 1,903 repeated-sequence rows were targeting guides. Rows sharing an identical sequence had different UIDs and different 12-sample count profiles. A count profile refers to the complete vector of experimental counts for one row.

Each Library B row was labelled using the `sequence_status` column:

- `unique_sequence`;
- `repeated_multi_gene`.

The resulting classification was:

| Guide type | Repeated multi-gene | Unique sequence |
| ---------- | ------------------: | --------------: |
| Non-targeting control | 0 | 1,000 |
| Targeting | 1,903 | 55,125 |

The repeated-sequence rows were retained because automatically excluding them would substantially reduce guide coverage:

- 362 genes would lose all guides;
- 277 genes would retain only one guide;
- 382 genes would retain two guides; and
- 1,021 genes would be affected in total.

Before sequence-based modelling, the origin of the repeated sequences and their distinct count profiles should be investigated. If the data are divided into training and testing sets, identical guide sequences must be assigned to the same partition to prevent information leakage.

## Processed Library B dataset

The complete cleaned and labelled Library B dataset was saved as:

```text
data/processed/library_b_prepared.rds
```
The processed file contains all 58,028 Library B rows. It includes:

- targeting guides and non-targeting controls;
- corrected `DEC1` annotations;
- validated 20-nucleotide guide sequences;
- all 12 experimental count columns;
- the `target_type` classification;
- sequence multiplicity information; and
- the `sequence_status` classification.

No control guides or repeated-sequence rows were permanently removed from this file.

The processed dataset can be loaded in R using:

```r
library_b_data <- readRDS(
  "data/processed/library_b_prepared.rds"
)
```
The targeting and control subsets should be regenerated from the latest version of `library_b_data`:

```r
library_b_targeting <- library_b_data[
  library_b_data$target_type == "targeting",
]

library_b_controls <- library_b_data[
  library_b_data$target_type == "non_targeting_control",
]
```


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

Open `Advanced-Bioinformatics.Rproj` in RStudio so that the repository root is used as the working directory.

Ensure that the original Excel files are available in:

```text
data/raw/
```

Install the required packages if they are not already installed:

```r
install.packages(c("readxl", "dplyr"))
```

Run the initial quality-control analysis:

```r
source("R/02_quality_control.R")
```

This script sources `R/01_import_data.R` automatically before performing the initial checks.

Prepare and save the complete Library B dataset:

```r
source("R/03_prepare_analysis_data.R")
```

This script:

1. imports the original data;
2. selects Library B;
3. joins counts and guide sequences by UID;
4. corrects the derived `DEC1` annotations;
5. validates sequences and counts;
6. classifies targeting guides and controls;
7. labels repeated guide sequences; and
8. saves the prepared dataset as:

```text
data/processed/library_b_prepared.rds
```

The saved dataset can be loaded in a later R session using:

```r
library_b_data <- readRDS(
  "data/processed/library_b_prepared.rds"
)
```

The analysis does not depend on objects stored in `.RData`.



