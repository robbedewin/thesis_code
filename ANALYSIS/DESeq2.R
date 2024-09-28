.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")

# Install Bioconductor and DESeq2 if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer"))

# Load the libraries
library(DESeq2)        # For differential expression analysis
library(ggplot2)       # For plotting
library(pheatmap)      # For heatmaps
library(RColorBrewer)  # For color palettes
library(dplyr)         # For data manipulation
library(tibble)        # For data frame manipulation
library(readr)         # For reading data
library(EnhancedVolcano) # For volcano plots (optional)
library(apeglm)        # For log fold change shrinkage




# Define file paths
counts_file <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/combined_counts_common_genes.tsv"
clinical_file <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/all_clinical_data.tsv"

# Load counts data
counts_data <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Load clinical data
clinical_data <- read.table(clinical_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert 'condition' and 'ETP_status' to factors
clinical_data$condition <- factor(clinical_data$condition)
clinical_data$ETP_status <- factor(clinical_data$ETP_status)


# Inspect the first few rows of counts data
head(counts_data)

# Inspect the first few rows of clinical data
head(clinical_data)

# Extract sample IDs from counts data
counts_samples <- colnames(counts_data)

# Extract sample IDs from clinical data
clinical_samples <- clinical_data$sample_id

# Check if all clinical samples are present in counts data
all_present <- all(clinical_samples %in% counts_samples)
print(paste("Are all clinical samples present in counts data?", all_present))  # Should return TRUE

# Identify any missing samples
missing_samples <- setdiff(clinical_samples, counts_samples)
print(paste("Number of missing samples:", length(missing_samples)))  # Should be 0

# If TRUE and no missing samples, proceed

# Reorder counts data columns to match the order of clinical data
counts_data <- counts_data[, clinical_samples]

# Verify alignment
all(colnames(counts_data) == clinical_data$sample_id)  # Should return TRUE

# Calculate row-wise sums to identify lowly expressed genes
gene_sums <- rowSums(counts_data)

# Plot distribution of gene counts
ggplot(data = data.frame(gene_sums), aes(x = gene_sums)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Distribution of Gene Counts", x = "Row-wise Counts (log10 scale)", y = "Frequency")

# save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/gene_counts_distribution.png")


# Define a threshold (e.g., keep genes with at least 10 counts in at least 20 samples)
keep_genes <- rowSums(counts_data >= 10) >= 20
sum(keep_genes)  # Number of genes retained

# Subset counts data
filtered_counts <- counts_data[keep_genes, ]

# Inspect dimensions
dim(filtered_counts)


# Create DESeq2 dataset for variance stabilizing transformation
dds_vst <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                  colData = clinical_data,
                                  design = ~ condition)  # Adjust the design formula based on your metadata

# Apply variance stabilizing transformation
vst_data <- vst(dds_vst, blind = FALSE)

# Perform PCA
pca <- plotPCA(vst_data, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

# Plot PCA
ggplot(pca, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of VST-transformed Counts") +
  theme_bw()

# save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/pca_plot.png")


# Calculate sample correlation matrix
sample_corr <- cor(assay(vst_data))

# Define color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot heatmap
pheatmap(sample_corr,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         col = colors,
         main = "Sample Correlation Heatmap")

# save the heatmap
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/sample_correlation_heatmap.png")



# Inspect unique conditions
unique(clinical_data$condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = clinical_data,
                              design = ~ condition)  # Adjust the design formula if you have more covariates

# Inspect the DESeq2 object
dds

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Obtain results for the condition comparison
# By default, DESeq2 compares the second level of the factor to the first
res <- results(dds)

# Inspect the results
summary(res)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# View the top differentially expressed genes
head(resOrdered)

# Shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = "condition_responsive_vs_IF", type = "apeglm")  # Adjust 'coef' based on your condition levels

# Inspect the shrunk results
head(resLFC)

# Save Volcano Plot
pdf("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/volcano_plot.pdf", width = 8, height = 6)
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                title = 'Volcano Plot',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0)
dev.off()


# Save MA Plot
pdf("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/ma_plot.pdf", width = 8, height = 6)
plotMA(resLFC, ylim = c(-5, 5), main = "MA Plot")
dev.off()


# Select top 20 DE genes
top_genes <- rownames(resOrdered)[1:20]

# Extract normalized counts
normalized_counts <- assay(vst(dds, blind = FALSE))[top_genes, ]

# Scale the data
scaled_counts <- t(scale(t(normalized_counts)))

# Define annotation for samples
annotation_col <- data.frame(
  Condition = clinical_data$condition
)
rownames(annotation_col) <- clinical_data$sample_id

# Save Heatmap
pdf("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/heatmap_top20_DE_genes.pdf", width = 10, height = 10)
pheatmap(scaled_counts,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         main = "Heatmap of Top 20 DE Genes")
dev.off()


# PCA plot using DESeq2's built-in function
plotPCA(vst_data, intgroup = c("condition")) + 
  theme_bw() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")

# Save PCA Plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/pca_plot.png", plot = last_plot(), width = 8, height = 6)


# Convert results to a data frame
res_df <- as.data.frame(resOrdered)

# Add gene symbols as a column
res_df$gene <- rownames(res_df)

# Save to TSV
write.table(res_df, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/results/deseq2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)




# Install and load clusterProfiler
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)  # Assuming human data

# Extract significant DE genes (adjusted p-value < 0.05 and |log2FC| > 1)
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(sig_genes, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

head(entrez_ids)
# Identify gene symbols with multiple Entrez IDs
duplicate_genes <- entrez_ids %>%
  group_by(SYMBOL) %>%
  filter(n() > 1) %>%
  pull(SYMBOL) %>%
  unique()

# Print the number of duplicated gene symbols
print(paste("Number of gene symbols with multiple Entrez IDs:", length(duplicate_genes)))

# View the duplicated gene symbols
duplicate_genes

# Option 2: Keep only the first Entrez ID per gene symbol
entrez_ids_unique <- entrez_ids %>%
  group_by(SYMBOL) %>%
  slice(1) %>%        # Keeps the first occurrence
  ungroup() %>%
  distinct(SYMBOL, ENTREZID)  # Ensures uniqueness

# Verify the uniqueness
any_duplicates_unique <- any(duplicated(entrez_ids_unique$SYMBOL))
print(paste("Are there any duplicate gene symbols after unique mapping?", any_duplicates_unique))  # Should return FALSE


# Extract unique Entrez IDs
entrez_ids_final <- entrez_ids_unique$ENTREZID

# Perform GO Enrichment Analysis
go_enrich <- enrichGO(gene = entrez_ids_final,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BF",          # BP: Biological Process MF: Molecular Function, CC: Cellular Component
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

# View GO enrichment results
head(go_enrich)

# Visualize GO enrichment with a dotplot
dotplot(go_enrich, showCategory = 20) + ggtitle("GO Enrichment of DE Genes")

# Save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/go_enrichment_dotplot_BF.png", plot = last_plot(), width = 10, height = 8)


# Perform KEGG enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID,
                          organism = 'hsa',  # 'hsa' for Homo sapiens
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# View enrichment results
head(kegg_enrich)

# Visualize KEGG enrichment
dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment of DE Genes")

# Save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/plots/KEGG_enrichment_dotplot.png", plot = last_plot(), width = 10, height = 8)


# Select significant DE genes
sig_de_genes <- rownames(res_df)[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)]

# Extract normalized counts for these genes
sig_normalized_counts <- assay(vst(dds, blind = FALSE))[sig_de_genes, ]

# Scale the data
scaled_sig_counts <- t(scale(t(sig_normalized_counts)))

# Define annotation for samples
annotation_col <- data.frame(
  Condition = clinical_data$condition
)
rownames(annotation_col) <- clinical_data$sample_id

# Plot heatmap
pdf("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/results/significant_DE_genes_heatmap.pdf", width = 10, height = 10)
pheatmap(scaled_sig_counts,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,  # Hide gene names for clarity
         show_colnames = FALSE,
         fontsize_row = 6,
         main = "Heatmap of Significant DE Genes")
dev.off()


