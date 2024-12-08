generate_plot <- function(airesultsfile, l2fcdf, filter_condition, output_dir) {
  # Load the allelic imbalance results
  airesults <- read.delim(file = airesultsfile, sep = "\t", as.is = TRUE)
  
  # Apply the filter condition
  airesults_sub <- subset(airesults, eval(parse(text = filter_condition)))
  
  # Subset the l2fcdf data
  l2fcdf_sub <- l2fcdf[l2fcdf$gene_name %in% airesults_sub$gene_name, ]
  l2fcdf_sub$gene_name <- factor(x = l2fcdf_sub$gene_name, levels = l2fcdf_sub[order(l2fcdf_sub$mean_expression, decreasing = FALSE), "gene_name"])
  
  # Load required libraries
  library(reshape2)
  library(ggplot2)
  
  # Create imbalanced gene-sample combinations
  imbalanced_gene_sample_combinations <- paste0(airesults_sub$gene_name, "_", airesults_sub$samples)
  
  # Melt the l2fcdf data
  l2fcdf_melt <- melt(data = l2fcdf_sub, id.vars = c("gene_name", "mean_expression"), variable.name = "sample_id", value.name = "l2fc")
  l2fcdf_melt$is_ai <- paste0(l2fcdf_melt$gene_name, "_", l2fcdf_melt$sample_id) %in% imbalanced_gene_sample_combinations
  
  # Create the plot
  p1 <- ggplot(data = l2fcdf_melt, mapping = aes(x = gene_name, y = l2fc))
  p1 <- p1 + geom_point(alpha = .5, shape = 16, size = 1.5, stroke = 0)
  p1 <- p1 + geom_violin(scale = "width", fill = "grey", alpha = .5)
  p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$is_ai, ], aes(x = gene_name, y = l2fc, fill = is_ai), shape = 21, size = 2, stroke = 0)
  p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  
  # Generate consistent file names based on the filter condition
  filter_name <- gsub(" ", "_", gsub(">", "GreaterThan", gsub("<", "SmallerThan", gsub("=", "EqualTo", filter_condition))))
  output_plot <- file.path(output_dir, paste0("l2fc_vst_", filter_name, ".png"))
  output_data <- file.path(output_dir, paste0("l2fc_vst_", filter_name, ".txt"))
  
  # Save the plot
  ggsave(filename = output_plot, plot = p1, width = 10, height = 5, dpi = 300)
  
  # Save the data
  write.table(x = l2fcdf_sub, file = output_data, quote = FALSE, sep = "\t", row.names = FALSE)
}

# Example usage:
generate_plot(
  airesultsfile = "/staging/leuven/stg_00096/home/rdewin/ASE/results/alloccurences_vst.txt",
  l2fcdf = l2fcdf,
  filter_condition = "n_up > 5",
  output_dir = "/staging/leuven/stg_00096/home/rdewin/ASE/results"
)
