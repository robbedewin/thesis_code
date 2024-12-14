.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

# Load the required libraries
library(DESeq2)
library(GenomicFeatures)
library(ggplot2)

# Load the functions
source(file = "/staging/leuven/stg_00096/home/rdewin/ASE/scripts/add_expression_data_functions.R")


# Make counts matrix for DESeq2, normalize and perform variance-stabilizing transformation
counts_file <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/counts/combined_counts.tsv"
matchedSamples <- c("P011", "P013", "P016", "P017", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P028", "P029", "P033", "P037", "P041", "P057", "P058", "P059", "P060", "P061", "P064", "P065", "P066", "P086", "P103", "P105")
output_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data"

l2fcfile <- process_expression_data(counts_file, matchedSamples, output_dir)
l2fcfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAlog2fc_vst.txt"
l2fcdf <- read.delim(file = l2fcfile, as.is = T)

## Combine p-values (of powered SNP loci) per gene and adjust for multiple testing

# Load gene annotation and log2-fold change data
hsexondb <- load_gene_annotations(gtffile = "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf")


for (SAMPLEID in matchedSamples) {
  print(paste("Processing sample:", SAMPLEID))
  
  # Overlap ASE results with gene annotations (exons) 
  ase_results_annot <- process_ase_results(SAMPLEID, l2fcdf, hsexondb) 
  # Combine p-values for same gene and adjust for multiple testing
  outdf <- combine_pvals_and_adjust(ase_results_annot)
  # Add log2 fold-change and gene names
  outdf <- add_log2fc_and_gene_names(outdf, l2fcdf, SAMPLEID)
  outdf <- format_output_data(outdf)
  save_output_data(outdf, output_dir, SAMPLEID)
  save_plot_imbalance(outdf, output_dir, SAMPLEID, significance = 0.01)
}


