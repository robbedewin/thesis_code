.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")

# Install Bioconductor and DESeq2 if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # For gene annotation


# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

combined_counts <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/counts/combined_counts.tsv"

# Import count data
count_data <- read.csv(combined_counts, sep = "\t", row.names = 1)

# Import sample metadata
col_data <- read.csv("path/to/sample_metadata.csv", row.names = 1)

# Ensure that sample names match between count data and metadata
all(colnames(count_data) == rownames(col_data))
# If FALSE, reorder metadata to match count data
col_data <- col_data[colnames(count_data), ]
