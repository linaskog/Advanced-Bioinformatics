# 04_sample_level_qc.R
# Advanced Bioinformatics: sample-level quality control
# Purpose:
# Evaluate sequencing depth, zero counts, sample correlations,
# normalized distributions and PCA for the prepared Library B data.

library(ggplot2)
library(DESeq2)
# Load the complete prepared Library B dataset
library_b_data <- readRDS(
  "data/processed/library_b_prepared.rds"
)

# Identify the 12 experimental count columns
count_columns <- grep(
  "^(BC-3|BJAB)_Day(0|14)_Rep[123]$",
  names(library_b_data),
  value = TRUE
)

count_columns

stopifnot(length(count_columns) == 12)

#total seq depth per sample
sequencing_depth<- data.frame(
  sample=count_columns,
  total_reads= colSums(library_b_data[, count_columns],
                       na.rm = TRUE),
  row.names = NULL
  
)
sequencing_depth

#QC figurs folder
dir.create(
  "results/qc",
  recursive = TRUE,
  showWarnings = FALSE
)


# seq_depth barplot
sequencing_depth$sample <- factor(
  sequencing_depth$sample,
  levels = count_columns
) 

sequencing_depth$group <- sub(
  "_Rep[123]$",
  "",
  as.character(sequencing_depth$sample)
)

sequencing_depth_plot<- ggplot(
  sequencing_depth, aes(
    x=sample, y=total_reads/1e6,
    fill=group
  )
)+
  geom_col(width = 0.75)+
  labs(
    title = "Library B sequencing depth by sample",
    x=NULL,
    y="Total reads(millions)",
    fill="Experimental group"
  )+
  theme_minimal(base_size = 12)+
  theme(
    axis.text.x = element_text(
      angle=45,
      hjust = 1
    )
  )
sequencing_depth_plot  

#save the plot
ggsave(
  filename = "results/qc/library_b_sequencing_depth.png",
  plot = sequencing_depth_plot,
  width = 10,
  height = 6,
  dpi = 300
)

file.exists("results/qc/library_b_sequencing_depth.png")

#zero_count percentage per sample
zero_counts_summary<- data.frame(
  sample= count_columns,
  zero_count_percentage= colMeans(library_b_data[, count_columns]==0, 
                                  na.rm = TRUE)*100,
row.names=NULL
)

zero_counts_summary


#zero count plot
zero_counts_summary$sample<- factor(
  zero_counts_summary$sample,
  levels = count_columns
)
zero_counts_summary$group<- sub(
  "_Rep[123]$",
  "",
  as.character(zero_counts_summary$sample)
)

zero_count_plot<- ggplot(
  zero_counts_summary, aes(
    x= sample,
    y= zero_count_percentage,
    fill = group)
)+
  geom_col(width = 0.75) +
  geom_text(
    aes(label = sprintf("%.2f", zero_count_percentage)),
    vjust = -0.4,
    size = 3
  ) +
  labs(
    title = "Zero-count guides in Library B samples",
    x = NULL,
    y = "Guides with zero counts (%)",
    fill = "Experimental group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

zero_count_plot

# Save the plot
ggsave(
  filename = "results/qc/library_b_zero_count_percentage.png",
  plot = zero_count_plot,
  width = 10,
  height = 6,
  dpi = 300
)

file.exists("results/qc/library_b_zero_count_percentage.png")

#correlation between replicates
count_matrix<- as.matrix(
  library_b_data[, count_columns]
)
log_count_matrix<- log2(count_matrix+1)
# Pearson correlations between samples
raw_correlation_matrix<- cor(
  log_count_matrix,
  method = "pearson")
round(raw_correlation_matrix,2)


# Convert the correlation matrix to long format for ggplot
correlation_plot_data<- as.data.frame(as.table(raw_correlation_matrix))

names(correlation_plot_data)<- c("sample_1",
"sample_2",
"correlation"
)

# Preserve the experimental sample order
correlation_plot_data$sample_1<- factor(
  correlation_plot_data$sample_1,
  levels = count_columns
)

correlation_plot_data$sample_2<- factor(
  correlation_plot_data$sample_2,
  levels = rev(count_columns)
)

#heat map
raw_correlation_heatmap<- ggplot(
  correlation_plot_data, aes(
    x= sample_1,
    y= sample_2,
    fill = correlation
  )
)+
  geom_tile(color="white")+
  geom_text(
    aes(label=sprintf("%.2f", correlation)),
    size= 2.7
  )+
  scale_fill_gradient(
    low="white",
    high = "darkblue",
    limits = c(0.5,1)
  )+
  coord_equal() +
  labs(
    title = "Library B sample correlations",
    subtitle = "Pearson correlation of log2(raw count + 1)",
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )
raw_correlation_heatmap

# Save the raw-count correlation heatmap
ggsave(
  filename = "results/qc/library_b_raw_correlation_heatmap.png",
  plot = raw_correlation_heatmap,
  width = 11,
  height = 9,
  dpi = 300
)

file.exists("results/qc/library_b_raw_correlation_heatmap.png")  




# Create sample metadata
sample_information <- data.frame(
  sample = count_columns,
  cell_line = ifelse(
    grepl("^BC-3", count_columns),
    "BC-3",
    "BJAB"
  ),
  day = ifelse(
    grepl("_Day0_", count_columns),
    "Day0",
    "Day14"
  ),
  replicate = sub(
    ".*_Rep",
    "Rep",
    count_columns
  ),
  row.names = count_columns
)

sample_information

#saveexperimental variable as factor
sample_information$cell_line <- factor(
  sample_information$cell_line,
  levels = c("BC-3", "BJAB")
)

sample_information$day <- factor(
  sample_information$day,
  levels = c("Day0", "Day14")
)

table(sample_information$day, useNA = "ifany")

rownames(count_matrix)<- library_b_data$UID
stopifnot(
  all(count_matrix >=0),
  all(count_matrix == floor(count_matrix))
)
storage.mode(count_matrix)<- "integer"

#deseq2 dataset
dds_qc<- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_information,
  design = ~ cell_line * day
)

#estimate sample specefic normalization factor
dds_qc<- estimateSizeFactors(dds_qc)
sizeFactors(dds_qc)


#variance-stabilizing transformation

vst_qc<- vst(dds_qc,
             blind = TRUE)
vst_matrix<- assay(vst_qc)
dim(vst_qc)

#correlation after normalizatiom and transformation
vst_matrix<-assay(vst_qc)
vst_correlation_matrix<- cor(
  vst_matrix,
  method = "pearson"
)

#BJAB day 14 replicate correlation
bjab_day14_samples<- grep(
  "^BJAB_Day14",
  count_columns,
  value = TRUE
)
round(
  vst_correlation_matrix[
    bjab_day14_samples,
    bjab_day14_samples],
  2
)


# heatmap after normalization
# Convert the VST correlation matrix to long format
vst_correlation_plot_data<- as.data.frame(
  as.table(vst_correlation_matrix)
)
names(vst_correlation_plot_data)<- c(
  "sample_1",
  "sample_2",
  "correlation"
)

vst_correlation_plot_data$sample_1<- factor(
  vst_correlation_plot_data$sample_1,
  levels = count_columns
)
vst_correlation_plot_data$sample_2<- factor(
  vst_correlation_plot_data$sample_2,
  levels = rev(count_columns)
)

vst_correlation_heatmap<- ggplot(
  vst_correlation_plot_data, aes(
    x=sample_1,
    y=sample_2,
    fill = correlation
  )
)+
  geom_tile(color="white")+
  geom_text(
    aes(label=sprintf("%.2f", correlation)),
    size = 2.7
  )+
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    limits = c(0.5,1)
  )+
  coord_equal()+
  labs(
    title="Library B normalized sample correlations",
    subtitle = "Pearson correlation after DESeq2 VST",
    x=NULL,
    y= NULL,
    fill="correlation"
  )+
  theme_minimal(base_size = 11)+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )
vst_correlation_heatmap

# Save the normalized-correlation heatmap
ggsave(
  filename = "results/qc/library_b_vst_correlation_heatmap.png",
  plot = vst_correlation_heatmap,
  width = 11,
  height = 9,
  dpi = 300
)

file.exists("results/qc/library_b_vst_correlation_heatmap.png")


#PCA 

pca_data <- plotPCA(
  vst_qc,
  intgroup= c("cell_line", "day"),
  ntop=500,
  returnData=TRUE
)
pca_data$replicate<- sample_information[
  rownames(pca_data),
  "replicate"
  ]
#the precentage 0f variance explained
pca_variance<- round(
  100*attr(pca_data, "percentVar"),
  1
)
pca_variance
pca_data


# Create clean experimental-group labels
pca_data$experimental_group <- factor(
  paste(pca_data$cell_line, pca_data$day),
  levels = c(
    "BC-3 Day0",
    "BC-3 Day14",
    "BJAB Day0",
    "BJAB Day14"
  )
)

# Create the PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = experimental_group,
    shape = replicate
  )
) +
  geom_point(size = 4) +
  labs(
    title = "PCA of Library B samples",
    subtitle = "DESeq2 VST using the 500 most variable guides",
    x = paste0("PC1: ", pca_variance[1], "% variance"),
    y = paste0("PC2: ", pca_variance[2], "% variance"),
    color = "Experimental group",
    shape = "Replicate"
  ) +
  theme_minimal(base_size = 12)

pca_plot


# Save the PCA plot
ggsave(
  filename = "results/qc/library_b_vst_pca.png",
  plot = pca_plot,
  width = 10,
  height = 7,
  dpi = 300
)

file.exists("results/qc/library_b_vst_pca.png")


# Convert log-transformed raw counts into long format
raw_distribution_data <- data.frame(
  sample = rep(
    count_columns,
    each = nrow(count_matrix)
  ),
  log2_count = as.vector(
    log2(count_matrix + 1)
  )
)

# Keep samples in experimental order
raw_distribution_data$sample <- factor(
  raw_distribution_data$sample,
  levels = count_columns
)

# Plot distributions before normalization
raw_distribution_plot <- ggplot(
  raw_distribution_data,
  aes(
    x = sample,
    y = log2_count,
    fill = sample
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  labs(
    title = "Library B count distributions before normalization",
    x = NULL,
    y = "log2(raw count + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

raw_distribution_plot
# Save the raw count-distribution plot
ggsave(
  filename = "results/qc/library_b_raw_count_distributions.png",
  plot = raw_distribution_plot,
  width = 11,
  height = 7,
  dpi = 300
)

file.exists("results/qc/library_b_raw_count_distributions.png")


# Convert VST-normalized values into long format
vst_distribution_data <- data.frame(
  sample = rep(
    count_columns,
    each = nrow(vst_matrix)
  ),
  vst_value = as.vector(vst_matrix)
)

# Keep samples in experimental order
vst_distribution_data$sample <- factor(
  vst_distribution_data$sample,
  levels = count_columns
)

# Plot distributions after normalization
vst_distribution_plot <- ggplot(
  vst_distribution_data,
  aes(
    x = sample,
    y = vst_value,
    fill = sample
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  labs(
    title = "Library B count distributions after normalization",
    subtitle = "DESeq2 variance-stabilizing transformation",
    x = NULL,
    y = "VST-transformed count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

vst_distribution_plot

# Save the VST-normalized distribution plot
ggsave(
  filename = "results/qc/library_b_vst_count_distributions.png",
  plot = vst_distribution_plot,
  width = 11,
  height = 7,
  dpi = 300
)

file.exists("results/qc/library_b_vst_count_distributions.png")

# Combine the main sample-level QC measurements
sample_qc_summary <- data.frame(
  sample = count_columns,
  total_reads = colSums(count_matrix),
  zero_count_percentage = colMeans(count_matrix == 0) * 100,
  deseq2_size_factor = sizeFactors(dds_qc)[count_columns],
  PC1 = pca_data[count_columns, "PC1"],
  PC2 = pca_data[count_columns, "PC2"],
  row.names = NULL
)

sample_qc_summary

# Save the sample-level QC summary
write.csv(
  sample_qc_summary,
  file = "results/qc/library_b_sample_qc_summary.csv",
  row.names = FALSE
)

file.exists("results/qc/library_b_sample_qc_summary.csv")
