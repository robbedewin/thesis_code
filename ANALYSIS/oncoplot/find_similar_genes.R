find_similar_genes <- function(target_genes, reference_genes, max_dist = 2) {
  similar_genes <- list()
  for (gene in target_genes) {
    distances <- stringdist::stringdist(gene, reference_genes)
    similar <- reference_genes[distances <= max_dist]
    if (length(similar) > 0) {
      similar_genes[[gene]] <- similar
    }
  }
  return(similar_genes)
}


# Extract gene names from genes_gr
reference_genes <- genes_gr$gene_name

# Add genes to genes_not_found (not deleting the previous ones)
genes_not_found <- c(genes_not_found, "H1F5", "NKX2.1")

# Find similar genes
similar_genes <- find_similar_genes(genes_not_found, reference_genes)

# Print the results
print(similar_genes)
