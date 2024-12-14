library(DESeq2)
library(ggplot2)

process_expression_data <- function(counts_file, matchedSamples, output_dir) {
  # Load counts matrix
  counts_data <- read.delim(counts_file, header = TRUE, row.names = 1, as.is = TRUE)
  
  # Filter counts dataframe to the matched samples
  filtered_counts_data <- counts_data[, colnames(counts_data) %in% matchedSamples]
  
  # Build a dummy sample table (all same condition, just to get normalized counts)
  sampleTable <- data.frame(sample = colnames(filtered_counts_data), condition = "T-ALL")
  
  # Create DESeq2 object from the filtered counts matrix
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts_data, 
                                colData   = sampleTable, 
                                design    = ~ condition)

  dds <- DESeq(dds)
  res <- results(dds)
  plotMA(res, ylim=c(-2, 2))
  # Estimate size factors (for normalization) 
  dds <- estimateSizeFactors(dds)
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)  # rows = genes, columns = samples
  
  #---------------------------------------------------------
  # STEP 1: Compute a "virtual mean sample" in log2 space
  #---------------------------------------------------------
  #   We add a small pseudocount (e.g., +1) to avoid log2(0).
  #   Then for each gene, average across all real samples to get a "virtual" sample's expression.
  log2_norm_counts <- log2(norm_counts + 1)
  virtual_mean_log2 <- rowMeans(log2_norm_counts)  # length = number of genes
  
  #---------------------------------------------------------
  # STEP 2: Compute log2 fold change for each sample vs virtual mean
  #---------------------------------------------------------
  #   log2FC_{sample} = log2(NormalizedCounts_sample + 1) - log2(VirtualMean + 1)
  #   That’s effectively: log2FC = log2_norm_counts - virtual_mean_log2
  log2fc_mat <- sweep(log2_norm_counts, 1, virtual_mean_log2, FUN = "-")
  
  #---------------------------------------------------------
  # STEP 3: Pick one sample to plot (analogous to left panel)
  #---------------------------------------------------------
  #   If you have multiple matchedSamples, pick one to illustrate; say the first column.
  sample_idx <- "P019"
  chosen_sample_name <- colnames(norm_counts)[sample_idx]
  
  # X-axis: the mean normalized count (on original scale or log scale).
  # Y-axis: log2 fold change (sample vs. virtual mean).
  df_plot <- data.frame(
    mean_norm_counts = rowMeans(norm_counts),           # average expression across all real samples
    log2fc           = log2fc_mat[, sample_idx]
  )
  
  # Optionally use log scale on x-axis (typical MA plot style).
  p0 <- ggplot(df_plot, aes(x = mean_norm_counts, y = log2fc)) +
    geom_point(alpha = 0.4) +
    scale_x_log10() +
    geom_hline(yintercept = 0, color = "red") +
    labs(
      title = paste("MA Plot: Sample", chosen_sample_name, "vs Virtual Mean"),
      x     = "Mean of Normalized Counts (log10 scale)",
      y     = "Log2 Fold Change vs Virtual Mean"
    )
  
  # Save the plot
  ggsave(filename = file.path(output_dir, "MAplot_virtual_mean_raw.png"), plot = p0, width = 6, height = 5)
  
  #---------------------------------------------------------
  # (OPTIONAL) Continue with VST and your prior code
  #---------------------------------------------------------
  dds_vst <- vst(dds, blind = TRUE)
  
  # Compute the "virtual mean sample" from the VST data
  dds_vstmeans <- rowMeans(assay(dds_vst))
  vst_fc <- assay(dds_vst) - dds_vstmeans
  
  # Plot VST-based log2FC vs mean expression
  p_vst <- ggplot(
      data = data.frame(
        mean_expression = dds_vstmeans, 
        log2fc = vst_fc[, sample_idx]
      ),
      aes(x = mean_expression, y = log2fc)
    ) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    labs(
      title = paste("VST-based MA: Sample", chosen_sample_name, "vs Virtual Mean"),
      x = "Mean Expression (VST)",
      y = "Log2 Fold Change (VST)"
    )
  
  ggsave(filename = file.path(output_dir, "MAplot_virtual_mean_vst.png"), plot = p_vst, width = 6, height = 5)
  
  #---------------------------------------------------------
  # Save tables if desired
  #---------------------------------------------------------
  # Write out the VST matrix
  write.table(assay(dds_vst), file.path(output_dir, "RNAcounts_vst.txt"), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  # Write out the raw normalized counts
  write.table(norm_counts, file.path(output_dir, "RNAcounts_normalised_T-ALL.txt"), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  # Write out the log2FC (VST) 
  l2fcdf <- as.data.frame(vst_fc)
  l2fcdf$gene_name <- rownames(l2fcdf)
  l2fcdf$mean_expression <- 2^rowMeans(log2(norm_counts + 1))
  write.table(l2fcdf, file.path(output_dir, "RNAlog2fc_vst.txt"), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  return(file.path(output_dir, "MAplot_virtual_mean_raw.png"))
}

# Example usage
process_expression_data(counts_file, matchedSamples, output_dir)
