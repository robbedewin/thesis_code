args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

.libPaths()
library(maftools)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(tidyverse)

# Define paths and load MAF files
maffiles <- list.files(path = input_dir, pattern = "\\.maf$", full.names = TRUE, recursive = TRUE)

# Read and merge MAF files
r_maffiles <- lapply(maffiles, read.maf)
unlistmaf <- maftools::merge_mafs(r_maffiles, verbose = TRUE)

# Generate plots
plot_file <- file.path(output_dir, "maftools_summary_plot.html")
pdf(file.path(output_dir, "oncoplot.pdf"))
oncoplot(maf = unlistmaf, top = 63, fontSize = 0.4, showTumorSampleBarcodes = TRUE)
dev.off()

# Save the plot to an HTML file
htmltools::save_html(plotmafSummary(unlistmaf, showBarcodes = TRUE), plot_file)

# Optionally print outputs or summary statistics
print(unlistmaf@summary)
