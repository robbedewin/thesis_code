.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

# Load the required libraries


## Check recurrence of allelically imbalanced + up/downregulated genes

# List the results files from the ASE analysis with added expression data (VST)
resultfiles <- list.files(
  path = "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data",
  pattern = "_imbalance_expression_vst.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Read in the results files
ai_results <- lapply(X = resultfiles, FUN = read.delim, as.is = TRUE)

# Combine the results
aidf <- do.call(rbind, ai_results)

# Add the sample name
aidf$sample <- rep(
  sub(
    pattern = "_imbalance_expression_vst.txt",
    replacement = "",
    x = basename(resultfiles)
  ),
  sapply(X = ai_results, nrow)
)

# Filter the results
aidf <- aidf[
  aidf$padj <= 0.05 &
    (aidf$log2fc >= 1 | aidf$log2fc <= -0.73),
]

# Filter out genes that are not of interest
aidf <- aidf[
  aidf$contig != "X" &
    !grepl(aidf$gene_name, pattern = "HLA-*") &
    !grepl(aidf$gene_name, pattern = "^LOC", perl = TRUE) &
    !grepl(aidf$gene_name, pattern = "^IG[HLK].*", perl = TRUE) &
    !grepl(aidf$gene_name, pattern = "^TR[ABDG][VCDJ].*", perl = TRUE),
]

imbalanced_gene_sample_combinations <- paste0(aidf$gene_name, "_", aidf$sample)

# Count the recurrence of genes
recurrent_genes <- sort(table(aidf$gene_name), decreasing = TRUE)
head(recurrent_genes, n = 25)

# Write the results to a file
write.table(
  x = aidf,
  file = "/staging/leuven/stg_00096/home/rdewin/ASE/results/allelic_imbalance_pooledsamples_vst.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
