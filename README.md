# Advanced Bioinformatics Project

## Project overview

This project analyses data from a genome-wide CRISPR-Cas9 knockout screen performed in two B-cell cancer cell lines:

- **BC-3:** a KSHV-positive primary effusion lymphoma cell line.
- **BJAB:** a KSHV-negative Burkitt lymphoma cell line.

The screen used the **GeCKO v2 sgRNA library**. The library is divided into Library A and Library B. Together, they generally contain six sgRNAs per protein-coding gene, with approximately three guides from each library.

The long-term objective is to develop a model that predicts sgRNA performance or quality from sgRNA sequence features. The precise response variable used to represent sgRNA quality has not yet been finalized.

## Biological principle

CRISPR-Cas9 is used to knock out genes in a population of cells. The abundance of each sgRNA is measured at the beginning of the experiment, Day 0, and after selection, Day 14.

If an sgRNA becomes less abundant at Day 14, disrupting its target may have reduced cell growth or survival. However, the observed change can be influenced by:

- sgRNA efficiency;
- biological importance of the target gene;
- cell-line-specific effects;
- sequencing depth;
- experimental noise.

Therefore, raw fold change cannot automatically be interpreted as sgRNA quality. The response variable and validation strategy must be selected carefully.

## Data source

The data originate from the supplementary material of Manzano et al. (2018):

- `41467_2018_5506_MOESM4_ESM.xlsx`
- `41467_2018_5506_MOESM6_ESM.xlsx`

Associated publication:

[Manzano et al., Nature Communications, 2018](https://doi.org/10.1038/s41467-018-05506-9)

The raw Excel files are stored locally in:

``` text
data/raw/
```

The raw data are excluded from Git using `.gitignore`. They must not be modified directly.

## Imported datasets

The script `R/01_import_data.R` imports four datasets.

| R object | Excel sheet | Contents | Dimensions |
|--------------|--------------|--------------------------------|-------------:|
| `gecko_counts` | `C. GECKO_BC3_BJAB` | sgRNA counts for BC-3 and BJAB | 123,411 × 14 |
| `gecko_library` | `F. GeCKOlibrary2` | sgRNA identifiers, sequences and gene annotations | 123,411 × 3 |
| `bc3_mageck` | `GECKO_BC-3` | Gene-level MAGeCK results for BC-3 | 21,915 × 14 |
| `bjab_mageck` | `GECKO_BJAB` | Gene-level MAGeCK results for BJAB | 21,915 × 14 |

The count table contains:

- BC-3 Day 0: three replicates;
- BC-3 Day 14: three replicates;
- BJAB Day 0: three replicates;
- BJAB Day 14: three replicates.

## Repository structure

``` text
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

1.  dataset dimensions;
2.  column data types;
3.  missing values;
4.  duplicated identifiers;
5.  agreement of identifiers between tables;
6.  agreement of gene annotations;
7.  negative and zero counts;
8.  sequencing depth;
9.  agreement between biological replicates.

### `03_prepare_analysis_data.R`

Prepares Library B for downstream analysis. The script:

1.  classifies guides as Library A or Library B from their UID prefixes;
2.  filters the count and sequence tables independently to Library B;
3.  joins the tables safely by UID;
4.  corrects the `DEC1` Excel-formatting inconsistency in derived annotation columns;
5.  validates guide sequences and experimental counts;
6.  classifies targeting guides and non-targeting controls;
7.  identifies and labels repeated guide sequences without removing them; and
8.  creates the complete cleaned Library B dataset used for subsequent QC and modelling preparation. \## Quality-control results

### Dimensions and data types

All four datasets were imported with the expected dimensions.

- UID, gene and guide-sequence columns were imported as character data.
- Experimental counts and MAGeCK statistics were imported as numeric data.
- BC-3 and BJAB MAGeCK tables had identical column structures.

### Missing values

No missing values were detected in any of the four imported datasets.

### Duplicate identifiers

No duplicated identifiers were found in:

- `gecko_counts$UID`;
- `gecko_library$UID`;
- `bc3_mageck$id`;
- `bjab_mageck$id`.

### UID coverage

Both sgRNA-level tables contained the same set of 123,411 UIDs:

- count UIDs missing from the sequence library: 0;
- library UIDs missing from the count table: 0.

Although the UIDs were present in both tables, their associated gene annotations were not always consistent.

## UID-to-gene mapping problem

Matching the count and sequence tables by UID initially produced:

| Mapping result            | Number of sgRNAs | Percentage |
|---------------------------|-----------------:|-----------:|
| Matching gene annotation  |           81,175 |     65.78% |
| Different gene annotation |           42,236 |     34.22% |

The mismatch was systematic rather than random.

### Library-specific results

| GeCKO library | Total guides | Matching guides | Mismatching guides | Match rate |
|--------------|-------------:|-------------:|--------------:|-------------:|
| Library A     |       65,383 |          23,150 |             42,233 |     35.41% |
| Library B     |       58,028 |          58,025 |                  3 |     99.99% |

Almost all UID-to-gene inconsistencies occurred in Library A.

The three apparent Library B mismatches involved the gene `DEC1`. Excel represented the same value differently:

- `1-Dec` in the count table;
- `42339` in the sequence-library table.

After correcting this formatting problem in derived variables, Library B had zero remaining UID-to-gene mismatches. The original raw data were not modified.

## Library selection and preparation

The UID prefix identifies the GeCKO sublibrary:

| UID prefix | Library   | Number of guides |
|------------|-----------|-----------------:|
| `HGLibA_`  | Library A |           65,383 |
| `HGLibB_`  | Library B |           58,028 |

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

| Guide type             | Number of guides |
|------------------------|-----------------:|
| Targeting guides       |           57,028 |
| Non-targeting controls |            1,000 |
| **Total**              |       **58,028** |

The targeting guides represent 19,049 annotated genes. The controls are identified by annotations beginning with:

``` text
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

| Guide type            | Repeated multi-gene | Unique sequence |
|-----------------------|--------------------:|----------------:|
| Non-targeting control |                   0 |           1,000 |
| Targeting             |               1,903 |          55,125 |

The repeated-sequence rows were retained because automatically excluding them would substantially reduce guide coverage:

- 362 genes would lose all guides;
- 277 genes would retain only one guide;
- 382 genes would retain two guides; and
- 1,021 genes would be affected in total.

Before sequence-based modelling, the origin of the repeated sequences and their distinct count profiles should be investigated. If the data are divided into training and testing sets, identical guide sequences must be assigned to the same partition to prevent information leakage.

## Processed Library B dataset

The complete cleaned and labelled Library B dataset was saved as:

``` text
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

``` r
library_b_data <- readRDS(
  "data/processed/library_b_prepared.rds"
)
```

The targeting and control subsets should be regenerated from the latest version of `library_b_data`:

``` r
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

- 42,236 sgRNAs would be removed;
- 81,175 sgRNAs would remain;
- 1,867 targets would disappear completely;
- the lost targets would include all 1,864 miRNA targets, plus `DEC1`, `PLN` and `SPHAR`.

The original publication did not report the miRNA-screen results because miRNA knockout analysis has additional biological complications.

The miRNA mismatches are part of the Library A mapping problem. They are not an independent source of mismatch.

## Count quality

No negative counts were detected in any experimental sample.

Zero counts were uncommon in most samples. BJAB Day 14 contained more zero-count guides:

| Experimental group | Percentage of zero counts |
|--------------------|--------------------------:|
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

$$\log_2(\text{count}+1)$$

The within-group replicate correlations were:

| Experimental group | Correlation range | Interpretation             |
|--------------------|------------------:|----------------------------|
| BC-3 Day 0         |         0.94–0.96 | Very strong agreement      |
| BC-3 Day 14        |         0.80–0.83 | Moderate-to-good agreement |
| BJAB Day 0         |         0.88–0.89 | Strong agreement           |
| BJAB Day 14        |         0.54–0.58 | Weak agreement             |

BJAB Day 14 showed the weakest replicate agreement and the largest proportion of zero counts. These samples will not be removed automatically, but their variability must be considered when defining the model response and evaluating model performance.

## Current conclusions

The following conclusions have been reached:

1.  The required datasets can be imported reproducibly.
2.  No missing values, duplicate identifiers or negative counts were detected.
3.  Sequencing depth differs between samples, so normalization is required.
4.  BC-3 and BJAB Day 0 replicates show strong agreement.
5.  BJAB Day 14 replicates show relatively weak agreement.
6.  Library B has a reliable UID-to-sequence mapping after correction of the `DEC1` Excel-formatting issue.
7.  Library A contains a systematic UID-to-gene mapping inconsistency and cannot yet be safely used for guide-level sequence modeling.
8.  The dataset includes 1,864 miRNA targets, which require separate consideration.
9.  No raw data have been modified or deleted.

## Genome mapping and annotation

To add genomic context to the CRISPR screen data, the targeting sgRNAs from Library B were mapped to the human reference genome hg19.

### Preparing sgRNA sequences

`R/06_prepare_genome_mapping.R` prepares the targeting guides for genome mapping.

The script:

- loads `data/processed/library_b_prepared.rds`;
- keeps only targeting guides;
- verifies that all guide sequences are 20 nt long and contain only A/C/G/T;
- creates a FASTA file using the guide UID as sequence identifier;
- writes the result to:

`data/processed/library_b_targeting_guides.fa`

The final FASTA contains 57,028 targeting sgRNAs.

### Genome indexing and Bowtie mapping

The hg19 reference genome was indexed with Bowtie 1.3.1 on HPC2N using SLURM.

The SLURM scripts are stored in:

- `slurm/build_hg19_index.slurm`
- `slurm/map_guides_to_hg19.slurm`

`build_hg19_index.slurm` builds the Bowtie index from `hg19.fa`.

`map_guides_to_hg19.slurm` maps all targeting sgRNAs against hg19 using exact matching only (`-v 0`) and reports all exact genomic matches (`-a`).

Mapping results:

- 57,028 guides processed
- 57,028 guides had at least one exact hg19 alignment
- 52,357 guides mapped uniquely
- 4,671 guides had multiple exact genomic matches
- 333,078 total exact alignments were reported

The uniquely mapped guide coordinates are stored in:

`data/processed/guide_hg19_unique_coordinates.tsv`

### Genome annotation

`R/07_annotate_genome_mapping.R` annotates the uniquely mapped guides using the UCSC hg19 RefGene annotation:

`data/reference/hg19.refGene.gtf.gz`

The annotation is imported with `rtracklayer` and represented as `GRanges` objects. Guide coordinates are also converted to `GRanges`, and `findOverlaps()` is used to compare sgRNA positions with transcript and exon regions.

Guides are classified as:

- **exon**: overlaps at least one annotated exon;
- **intron**: overlaps an annotated transcript but no exon;
- **intergenic**: overlaps no annotated transcript.

Annotation summary:

- 52,281 exon-overlapping guides
- 26 intron-only guides
- 50 intergenic guides

Gene annotation is also added. If a guide overlaps more than one annotated gene, all distinct gene names are retained as a semicolon-separated value rather than choosing one arbitrarily.

Gene-overlap summary:

- 48,377 guides overlap exactly one gene
- 3,930 guides overlap multiple genes
- 50 guides overlap no annotated gene

The final annotated table is stored in:

`data/processed/guide_hg19_annotated.tsv`

The table contains genomic coordinates, strand, mapped sequence, region type, gene annotation, and the number of overlapping genes for each uniquely mapped sgRNA.

### BC-3 H3K27ac ChIP-seq integration

As a second additional dataset, BC-3 H3K27ac ChIP-seq peak data were integrated to add BC-3-specific regulatory chromatin information to the CRISPR guide annotations.

The processed peak file used was:

`data/raw/GSM4040945_BC3.merge.peak.bed.gz`

This BED file contains 25,546 H3K27ac-enriched genomic regions from the BC-3 cell line.

`R/08_integrate_bc3_h3k27ac.R`:

- imports the H3K27ac peak BED file with `rtracklayer`;
- converts the peak coordinates to a `GRanges` object;
- loads the uniquely mapped and genome-annotated sgRNA table;
- converts guide coordinates to `GRanges`;
- uses `findOverlaps()` to identify sgRNAs whose genomic target sites overlap BC-3 H3K27ac peaks;
- adds a binary feature, `H3K27ac_peak_overlap`.

Integration summary:

- 52,357 uniquely mapped guides analysed
- 25,546 BC-3 H3K27ac peaks
- 4,118 guides overlap at least one H3K27ac peak
- 48,239 guides do not overlap an H3K27ac peak

The integrated output is stored in:

`data/processed/guide_hg19_bc3_h3k27ac.tsv`

## The resulting table combines genomic coordinates, gene and region annotation, and BC-3 H3K27ac regulatory information for downstream statistical or machine-learning analyses.

## Reproducing the current analysis

Open `Advanced-Bioinformatics.Rproj` in RStudio so that the repository root is used as the working directory.

Ensure that the original Excel files are available in:

``` text
data/raw/
```

Install the required packages if they are not already installed:

``` r
install.packages(c("readxl", "dplyr"))
```

Run the initial quality-control analysis:

``` r
source("R/02_quality_control.R")
```

This script sources `R/01_import_data.R` automatically before performing the initial checks.

Prepare and save the complete Library B dataset:

``` r
source("R/03_prepare_analysis_data.R")
```

This script:

1.  imports the original data;
2.  selects Library B;
3.  joins counts and guide sequences by UID;
4.  corrects the derived `DEC1` annotations;
5.  validates sequences and counts;
6.  classifies targeting guides and controls;
7.  labels repeated guide sequences; and
8.  saves the prepared dataset as:

``` text
data/processed/library_b_prepared.rds
```

The saved dataset can be loaded in a later R session using:

``` r
library_b_data <- readRDS(
  "data/processed/library_b_prepared.rds"
)
```

The analysis does not depend on objects stored in `.RData`.

## Data Preparation

Following the quality control and annotation, the GeCKO library b sgRNA's were loaded along with the genome annotation data. The data points corresponding to non-targeting control sgRNAs were removed along with multimapping sgRNAs, resulting in 51907 data points for modelling.

The log2 Fold Change between Day 0 and Day 14 BC-3 cell line, sgRNA abundance counts was computed using the DESeq2 library. The gene-level average of log2 Fold Change was subtracted from these values to account for differences in sgRNA abundance between genes. These values were then used as the dependent variable to estimate sgRNA efficacy.

To model the sgRNA efficacy, the guide sequence composition was used. The overall GC-content was computed as the ratio of Guanine (G) and Cytosine (C) in the guide. Position-wise nucleotide information was also used as input for the model.

From the annotation experiment, the strand, chromosome, region type, gene overlap count and H3K27ac peak overlap information was extracted.

The variables were saved for downstream modelling as:

```         
data/processed/variables.rds
```

The saved dataset can be accessed with:

``` R
variables <- readRDS(
  "data/processed/variables.rds"
)
```

This data was subseqeuntly divided into training data (80%) and testing data (20%), separated with respect to gene names to ensure all sgRNAs for the same gene end up in the same set.

The data preprocessing procedure can be found in: "R/09_data_preparation.R"

## Linear Model

A base R linear regression model was trained to predict centered log2 fold change in sgRNA abundance using the previously prepared variables. This model

### Training

A linear model was trained to predict the centered log2 fold change in sgRNA abundance based on position-wise nucleotide data, GC-content, strand, chromosome, region type, gene overlap count and H3k27ac peak overlap data. The first model produced the following result:

``` text
Residual standard error: 0.3255 on 41459 degrees of freedom
Multiple R-squared:  0.04302,	Adjusted R-squared:  0.04152 
F-statistic: 28.67 on 65 and 41459 DF,  p-value: < 2.2e-16
```

The first model only explained 4.3% of the variance in the training data set. As a previous study optimizing sgRNA design found position-wise dinucleotide data contributes to activity predictions (Doench et al. 2016), we added cross-terms of the position-wise nucleotide variables in an attempt to improve model performance. The following result was obtained for the new model:

``` text
Residual standard error: 0.3236 on 41288 degrees of freedom
Multiple R-squared:  0.05838,	Adjusted R-squared:  0.053 
F-statistic: 10.85 on 236 and 41288 DF,  p-value: < 2.2e-16
```

The new model could explain 5.8% of the variation in the training data.

### Exploration

|              |              |                |             |                 |
|--------------|--------------|----------------|-------------|-----------------|
| **Variable** | **Estimate** | **Std. Error** | **t value** | **Pr(\>\|t\|)** |
| pos19_G      | -0.18580334  | 0.013325345    | -13.943605  | 4.333181e-44    |
| pos20_G      | -0.13073050  | 0.009563967    | -13.669067  | 1.922580e-42    |
| pos19_C      | -0.10780585  | 0.012595966    | -8.558760   | 1.179549e-17    |
| pos18_C      | -0.09885451  | 0.011576048    | -8.539573   | 1.392436e-17    |
| pos17_G      | -0.10235272  | 0.012542660    | -8.160367   | 3.433615e-16    |
| pos17_C      | -0.09508607  | 0.012221485    | -7.780238   | 7.406420e-15    |
| pos18_G      | -0.08948371  | 0.012638102    | -7.080470   | 1.459552e-12    |
| pos18_A      | -0.08231528  | 0.012703086    | -6.479944   | 9.278779e-11    |
| pos19_A      | -0.09102501  | 0.014078891    | -6.465354   | 1.021876e-10    |
| pos17_A      | -0.06792672  | 0.013530278    | -5.020349   | 5.179183e-07    |

## Extreme Gradient Boosting Model

XGBoost (Extreme Gradient Boosting) was used to model centered log2 fold change of sgRNA counts by sequence composition and genome annotation information. XGBoost is a scalable gradient boosting framework. The model builds sequential decision trees, with each tree attempting to improve the predictions of the previous trees.

The response variable is the centered log2 fold change in sgRNA abundance and the predictor variables the generated guide sequence and annotation information.

The training data was used to fit the model and for cross-validation during hyperparameter optimizations, whereas the testing data was only used for the evaluation of the final model.

### Hyperparameter Selection

Hyperparameters eta, gamma, max_depth, subsample, colsample_bytree, and min_child_weight were optimized for a small subset of values to improve model performance. Hyperparameter optimization was performed using 5-fold cross-validation on the training data. The number of boosting rounds was also determined during cross-validation using early stopping.

For each parameter configuration up to 1,000 boosting rounds were evaluated, stopping training when RMSE failed to improve for 30 consecutive rounds. The parameters were optimized sequentially in the previously stated order, retaining the selected parameter value from the previous optimization.

The cross-validation Root Mean Squared Error (CV RMSE) was used to select parameter values, where a lower value indicates a lower prediction error.

#### eta

eta or learning rate sets the step size shrinkage used in update to prevent overfitting. The default value for eta for a Tree Booster model is 0.3 and the parameter has a range of 0 to 1. For the optimization three values larger than 0.3 (0.4, 0.5 and 0.6) and four values smaller than 0.3 (0.15, 0.1, 0.05 and 0.01) were investigated.

Reported for each tested value is the cross-validation RMSE and Mean Absolute Error (MAE).

|         |            |            |           |
|---------|------------|------------|-----------|
| **eta** | **rounds** | **cvRMSE** | **cvMAE** |
| 0.01    | 254        | 0.3254402  | 0.2372527 |
| 0.05    | 78         | 0.3252503  | 0.2374313 |
| 0.10    | 40         | 0.3254807  | 0.2376654 |
| 0.15    | 25         | 0.3255142  | 0.2377236 |
| 0.30    | 6          | 0.3262934  | 0.2378328 |
| 0.40    | 4          | 0.3264669  | 0.2379760 |
| 0.50    | 3          | 0.3265922  | 0.2380111 |
| 0.60    | 2          | 0.3268551  | 0.2383120 |

Out of the tested values, eta = 0.05 produced the lowest CV RMSE and was used for subsequent optimizations.

#### gamma

gamma or min_split_loss refers to the minimum loss reduction required to make a split. A larger value makes the model more conservative. The default value for a Tree Booster is 0 and the range of the parameter is 0 to infinity. Gamma values 0, 0.05, 0.1, 0.15, 0.2 and 0.3 were investigated.

|           |            |            |           |
|-----------|------------|------------|-----------|
| **gamma** | **rounds** | **cvRMSE** | **cvMAE** |
| 0.00      | 78         | 0.3252503  | 0.2374313 |
| 0.05      | 73         | 0.3251241  | 0.2372781 |
| 0.10      | 73         | 0.3252413  | 0.2373863 |
| 0.15      | 81         | 0.3251851  | 0.2373946 |
| 0.20      | 75         | 0.3252455  | 0.2373683 |
| 0.30      | 71         | 0.3252376  | 0.2373063 |

Out of the tested values, gamma = 0.05, generated the lowest CV RMSE and was selected for subsequent optimizations.

#### max_depth

The max_depth hyperparameter controls the maximum depth of a tree. Increasing this value causes the model to become more complex and more likely to overfit the data. The default value for max_depth for a Tree Booster is 6, with the range 0 to infinity. Investigated values of max_depth were 3, 4, 5, 6, 7, 8 and 9.

|               |            |            |           |
|---------------|------------|------------|-----------|
| **max_depth** | **rounds** | **cvRMSE** | **cvMAE** |
| 3             | 230        | 0.3247622  | 0.2370123 |
| 4             | 147        | 0.3248728  | 0.2370525 |
| 5             | 103        | 0.3249376  | 0.2371689 |
| 6             | 73         | 0.3251241  | 0.2372781 |
| 7             | 60         | 0.3258134  | 0.2378620 |
| 8             | 38         | 0.3264278  | 0.2381419 |
| 9             | 25         | 0.3274354  | 0.2386453 |

From the tested values, max_depth = 3 generated the lowest CV RMSE and was selected for subsequent optimizations.

#### subsample

This parameter controls the subsample ratio of the training data. A value of 0.5 means that the model would randomly sample half the data prior to growing trees. The default value for a Tree Booster is subsample = 1 with a range of 0 to 1. Selected for investigation were values 1, 0.9, 0.8, 0.7, 0.6 and 0.5.

|               |            |            |           |
|---------------|------------|------------|-----------|
| **subsample** | **rounds** | **cvRMSE** | **cvMAE** |
| 1.0           | 230        | 0.3247622  | 0.2370123 |
| 0.9           | 214        | 0.3245957  | 0.2369308 |
| 0.8           | 203        | 0.3246269  | 0.2369205 |
| 0.7           | 194        | 0.3246410  | 0.2369731 |
| 0.6           | 172        | 0.3247639  | 0.2369658 |
| 0.5           | 164        | 0.3247871  | 0.2370838 |

From the tested values, subsample = 0.9 generated the lowest CV RMSE and was selected for subsequent optimizations.

#### colsample_bytree

The colsample_bytree parameter is the subsample ratio of columns when constructing each tree. The default value for a Tree Booster is 1 with a range of 0 to 1. Selected for investigation were values 1, 0.9, 0.8, 0.7, 0.6 and 0.5.

|                      |            |            |           |
|----------------------|------------|------------|-----------|
| **colsample_bytree** | **rounds** | **cvRMSE** | **cvMAE** |
| 1.0                  | 214        | 0.3245957  | 0.2369308 |
| 0.9                  | 211        | 0.3245992  | 0.2369133 |
| 0.8                  | 251        | 0.3245433  | 0.2369460 |
| 0.7                  | 210        | 0.3246796  | 0.2369299 |
| 0.6                  | 237        | 0.3245733  | 0.2369516 |
| 0.5                  | 236        | 0.3246678  | 0.236973  |

From the tested values, colsample_bytree = 0.8 generated the lowest CV RMSE and was selected for subsequent optimizations.

#### min_child_weight

This parameter determines the minimum sum of instance weight needed in a child. For a Tree Booster the default value i 1 with a range of 0 to infinity, where 0 indicates no limit. Selected for investigation were values 0, 0.5, 1, 2, 5 and 10.

|                      |            |            |           |
|----------------------|------------|------------|-----------|
| **min_child_weight** | **rounds** | **cvRMSE** | **cvMAE** |
| 0.0                  | 251        | 0.3245433  | 0.2369460 |
| 0.5                  | 251        | 0.3245433  | 0.2369460 |
| 1.0                  | 251        | 0.3245433  | 0.2369460 |
| 2.0                  | 195        | 0.3246925  | 0.2369244 |
| 5.0                  | 202        | 0.3246511  | 0.2369114 |
| 10.0                 | 209        | 0.3246218  | 0.2368987 |

From the tested values, min_child_weight 0, 0.5 and 1 generated identical CV RMSE values.

#### Selected Model Hyperparameters

The sequential optimization resulted in the following selected parameters:

| **Parameter**    | **Selected Value** |
|------------------|--------------------|
| eta              | 0.05               |
| gamma            | 0.05               |
| max_depth        | 3                  |
| subsample        | 0.9                |
| colsample_bytree | 0.8                |
| min_child_weight | 0                  |

After selecting these parameter values, a final cross-validation identified 251 boosting rounds.

The optimization procedure can be found in: "R/12_xgb_model_optimization.R"

### Final Model Evaluation

The final XGBoost model was trained on the complete training dataset using the selected hyperparameters and 251 boosting rounds. The held-out test dataset was then used to obtain an independent evaluation of model performance. The package basicstats was used to compute RMSE, MAE, R2 and Spearman's correlation coefficients between the observed and predicted datasets for the training data and the held-out testing data. The result of this analysis is shown in the table below.

|  |  |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|----|
| **rounds** | **trainRMSE** | **testRMSE** | **trainMAE** | **testMAE** | **trainR2** | **testR2** | **trainSpearman** | **testSpearman** |
| 251 | 0.3216164 | 0.3263604 | 0.234871 | 0.2370772 | 0.06428116 | 0.04966186 | 0.2423636 | 0.2198846 |

Measuring RMSE between observed and predicted values, the model performed similarly for the training and testing data, 0.32 and 0.33 respectively. The calculated MAE values were also similar between the datasets, 0.23 and 0.24 for the training and testing data respectively. The calculated Spearman correlation between observed and predicted values was 0.24 for the training-set and 0.22 for the testing-set. The computed R-squared for the final model was 0.05 for the test-data, whereas the R-squared for the training-data was 0.06.

The model showed positive association between predicted and observed sgRNA efficacy, and prediction errors were similar for both the training and testing data. The code for the final model evaluation can be found in: "R/13_xgb_model_evaluation.R"

## SHAP

![](results/models/shap_plot_xgb_model.pdf)

## Final Model Exploration

Differences in model predictions using the linear and XGBoost model were visualized in scatterplots, using the training and testing data. Predictions are shown to have an approximate 1-to-1 relationship. Furthermore, correlation between the linear model and XGBoost model predictions was \~0.9 for both the training and test data sets. And correlation between the absolute residuals between the predictions and observed values were 0.99.
