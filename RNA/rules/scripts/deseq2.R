.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")
library(DESeq2)
library(scales)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("/lustre1/project/stg_00096/home/rdewin/RNA")

#variables for testing outside of snakemake
input_counts <- "results/counts/combined_counts.tsv"
sample_info_file <- "results/sample_info.csv"
results_dds <- "results/deseq2/deseq2_results.tsv"

# Snakemake variables
#Input
input_counts <- snakemake@input$counts
sample_info_file <- snakemake@input$sample_info_file

#Output
results_dds <- snakemake@output[[1]]
output_dir <- dirname(results_dds)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#Params
plot_dir <- file.path(output_dir, "plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


# Open the log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

# Read in the count data and sample information
count_data <- read.csv(input_counts, sep="\t", row.names=1)
sample_info <- read.csv(sample_info_file, row.names = 1)
sample_info$ETPstatus <- factor(sample_info$ETPstatus)

# Create DESeq2 dataset with ETPstatus as the design
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ ETPstatus)

dds <- dds[rowSums(counts(dds)) > 10,] # Filter low count genes
dds <- DESeq(dds)  # Run DESeq analysis
results <- results(dds)

# Save DESeq2 results
write.csv(as.data.frame(results), file = results_dds)

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

# PCA plot with explained variance
pcaData <- plotPCA(vsd, intgroup="ETPstatus", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=ETPstatus)) +
    geom_point(size=3) +
    geom_text_repel(aes(label=rownames(pcaData))) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance"))
ggsave(file.path(plot_dir, "PCA_plot.png"), plot=p)

# Heatmap of sample-to-sample distances
# Visualize results with a heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$ETPstatus
colnames(sampleDistMatrix) <- vsd$ETPstatus
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         filename=file.path(plot_dir, "heatmap.png"))

# Volcano plot of results
volcanoData <- as.data.frame(results)
volcanoData$log2FoldChange <- -log10(volcanoData$pvalue) * sign(volcanoData$log2FoldChange)
ggvolcano <- ggplot(volcanoData, aes(x=log2FoldChange, y=-log10(pvalue), color=padj < 0.05)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=c("red", "grey")) +
    theme_minimal()
ggsave(file.path(plot_dir, "volcano_plot.png"), plot=ggvolcano)

# Convert gene IDs if necessary
sig_genes <- rownames(subset(results, padj < 0.05))
sig_genes_entrez <- mapIds(org.Hs.eg.db, keys=sig_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")

# GO enrichment analysis
ego <- enrichGO(gene = sig_genes_entrez,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                readable = TRUE)

# KEGG pathway analysis
ekegg <- enrichKEGG(gene = sig_genes_entrez,
                    organism = 'hsa')

# Visualize results
dotplot(ego)
dotplot(ekegg)

# Save plots
#ggsave(file.path(plot_dir, "ego_dotplot.png"), plot=dotplot(ego))
#ggsave(file.path(plot_dir, "ekegg_dotplot.png"), plot=dotplot(ekegg))
