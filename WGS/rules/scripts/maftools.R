.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")
library(maftools)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(tidyverse)

setwd("/staging/leuven/stg_00096/home/rdewin/zzz")
input_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/mutect2/"
output_html <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/maftools/maftools_summary.html"
output_pdf <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/maftools/maftools_oncoplot2.pdf"
output_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/maftools/"
# Define paths directly from Snakemake variables (set these in the Snakemake rule)
# input_dir <- snakemake@params[["dir"]]
# output_dir <- dirname(snakemake@output[["html"]])
# output_html <- snakemake@output[["html"]]
# output_pdf <- snakemake@output[["pdf"]]

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Open the log file
# log <- file(snakemake@log[[1]], open = "wt")
# sink(log)
# sink(log, type = "message")




# List the MAF files
maffiles <- list.files(path = input_dir, pattern = "P.*_annotated_lifted_variants\\.maf$", full.names = TRUE, recursive = TRUE)
if (length(maffiles) == 0) stop("No MAF files found in the directory")


# Read MAF files and merge them
r_maffiles <- lapply(maffiles, read.maf)
merged_maf <- maftools::merge_mafs(r_maffiles, verbose = TRUE)

# Read MAF files and merge them, assigning sample names based on file names
r_maffiles <- lapply(maffiles, function(maf_file) {
  maf_data <- read.maf(maf_file)
  sample_name <- tools::file_path_sans_ext(basename(maf_file))
  sample_name <- sub("^(.*?)_.*", "\\1", sample_name)
  maf_data@data$Tumor_Sample_Barcode <- sample_name  # Ensuring unique barcodes
  return(maf_data)
})
merged_maf <- maftools::merge_mafs(r_maffiles, verbose = TRUE)

# Print summary of the merged MAF object
print(merged_maf@summary)

# Generate Oncoplot
pdf(file.path(output_pdf))
oncoplot(maf = merged_maf, top = 25, fontSize = 0.8,
         showTumorSampleBarcodes = TRUE)
dev.off()

library(htmltools)
# Generate HTML summary report
maf_summary_plot <- plotmafSummary(maf = merged_maf, rmOutlier = TRUE, addStat = 'median', dashboard = FALSE)
maf_summary_plot

# Save the plot to an HTML file
save_html(maf_summary_plot, file = output_html)

# Optionally print outputs or summary statistics
print(merged_maf@summary)

# Close logging
sink()
sink(type = "message")