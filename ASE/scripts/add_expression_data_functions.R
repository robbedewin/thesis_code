## Functions
combine_pvals <- function(ase_pergene) {
  # ase_pergene <- ase_results_annot[4:14, ]
  outdf <- data.frame(contig = ase_pergene[1, "contig"], positions = paste0(unique(ase_pergene$position), collapse = ","), pcombined = 1, gene = ase_pergene[1, "gene"], stringsAsFactors = F, power = T)
  
  is_duplicated <- duplicated(ase_pergene$position)
  has_power <- ase_pergene$filter <= 0.02
  ase_pergene <- ase_pergene[!is_duplicated & has_power, ]
  
  if (nrow(ase_pergene) > 1) {
    outdf$pcombined <- fishersMethod(x = ase_pergene$pval)
  } else if (nrow(ase_pergene) == 1) {
    outdf$pcombined <- ase_pergene$pval
  } else {
    outdf$power <- F
  }
  return(outdf)
}

fishersMethod <- function(x) {
  pchisq(q = -2 * sum(log(x)), df = 2*length(x), lower.tail = F)
}

# Plot imbalance expression, with colouring for upregulated and downregulated genes
plot_imbalance_expression <- function(imbalancedf, significance = 0.01) {
  imbalancedf <-  imbalancedf[order(imbalancedf$mean_expression, decreasing = F), ]
  labeldf <- data.frame(pos = unlist(lapply(X = 10^(0:4), FUN = function(x) sum(imbalancedf$mean_expression < x))), expr = 10^(0:4), stringsAsFactors = F)
  
  outdf_bak <- imbalancedf
  imbalancedf <- imbalancedf[!grepl(pattern = "^HLA.*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^IG[HLK].*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^LOC", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^TR[ABDG][VCDJ].*", x = imbalancedf$gene_name, perl = T), ]

  imbalancedf$notes <- ifelse(imbalancedf$padj > significance, "nonsig", 
                              ifelse(imbalancedf$log2fc >= 1, "up",
                                     ifelse(imbalancedf$log2fc <= -.73, "down", "nonsig")))
  
  p1 <- ggplot(data = imbalancedf, mapping = aes(x = 1:nrow(imbalancedf), y = -sign(log2fc)*log10(pcombined)))
  p1 <- p1 + geom_point(mapping = aes(colour = notes, size = abs(log2fc)), 
                        show.legend = F, alpha = .4)
  p1 <- p1 + geom_hline(yintercept = c(-1,1)*-log10(max(imbalancedf[imbalancedf$padj < significance, "pcombined"])), linetype = "dashed", colour = "grey") +
    geom_text(data = imbalancedf[imbalancedf$notes != "nonsig", ], mapping = aes(x = which(imbalancedf$notes != "nonsig"), y = -sign(log2fc)*log10(pcombined), label = gene_name), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(imbalancedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_y_continuous(breaks = seq(-10,10,2), oob = scales::squish, limits = c(-10,10))
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$expr, name = "mean expression (normalised)")
  # p1 <- p1 + scale_color_brewer(type = "div", palette = "RdBu", direction = -1)
  p1 <- p1 + scale_color_manual(values = c(nonsig = "#e0e0e0", up = "#ef8a62", down = "#67a9cf"))
  p1 <- p1 + scale_size_continuous(range = c(1,7.5))
  p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = -90)) + labs(x = NULL)
  p1 <- p1 + annotate("text", x = -Inf, y = -Inf, label = paste("p-value =", significance), hjust = -0.1, vjust = -1.5, size = 3, color = "black", alpha = 0.5)
  
  return(p1)  
}


process_expression_data <- function(counts_file, matchedSamples, output_dir) {
  # Load counts matrix
  counts_data <- read.delim(counts_file, header = TRUE, row.names = 1, as.is = TRUE)
  
  # Filter counts dataframe
  filtered_counts_data <- counts_data[, colnames(counts_data) %in% matchedSamples]
  
  sampleTable <- data.frame(sample = colnames(filtered_counts_data), condition = "T-ALL")
  
  # Create DESeq2 object from counts matrix
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts_data, 
                                colData = sampleTable, 
                                design = ~1)
  
  
  # Normalize and perform variance-stabilizing transformation (VST)
  dds <- estimateSizeFactors(dds)
  dds_vst <- vst(dds, blind = TRUE)
  
  # Compute mean expression and log2 fold-change relative to "virtual mean sample"
  dds_vstmeans <- rowMeans(assay(dds_vst))
  vst_fc <- assay(dds_vst) - dds_vstmeans
  
  # Create a ggplot object for mean expression vs log2 fold-change
  p1 <- ggplot(data = data.frame(mean_expression = dds_vstmeans, log2fc = vst_fc[, 1]), 
         aes(x = mean_expression, y = log2fc)) +
    geom_point() +
    labs(title = "Log2 Fold Change vs Mean Expression", x = "Mean Expression", y = "Log2 Fold Change")

  # Save the plot
  ggsave(filename = file.path(output_dir, "mean_expression_vs_log2fc.png"), plot = p1)

  # Create a ggplot object for log10 mean normalized counts vs log2 fold-change
  p2 <- ggplot(data = data.frame(mean_normalized_counts = log10(rowMeans(counts(dds, normalized = TRUE)) + 1), 
                   log2fc = vst_fc[, 1]), 
         aes(x = mean_normalized_counts, y = log2fc)) +
    geom_point() +
    labs(title = "Log2 Fold Change vs Log10 Mean Normalized Counts", x = "Log10 Mean Normalized Counts", y = "Log2 Fold Change")

  # Save the plot
  ggsave(filename = file.path(output_dir, "log10_mean_normalized_counts_vs_log2fc.png"), plot = p2)
  
  # Save normalized and VST-transformed data
  write.table(assay(dds_vst), file.path(output_dir, "RNAcounts_vst.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  # Save log2 fold-change data
  l2fcdf <- as.data.frame(vst_fc)
  
  # Create normalized counts dataframe
  resdf <- as.data.frame(counts(dds, normalized = TRUE))
  
  # Add gene names (genes were already annotated in the counts file) and mean expression
  l2fcdf$gene_name <- rownames(l2fcdf)
  l2fcdf$mean_expression <- 2^rowMeans(log2(resdf + 1))

  # Save log2 fold-change data
  l2fcfile <- file.path(output_dir, "RNAlog2fc_vst.txt")
  write.table(l2fcdf, file = l2fcfile, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  # Save normalized counts
  write.table(resdf, file.path(output_dir, "RNAcounts_normalised_T-ALL.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

  return(l2fcfile)
}


load_gene_annotations <- function(gtffile) {
  library(GenomicFeatures)
  hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
  seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb))) # converts chr1 to 1
  hsexondb <- exons(x = hstxdb, columns = c("gene_id"))
  return(hsexondb)
}

process_ase_results <- function(SAMPLEID,
                                l2fcdf,
                                hsexondb,
                                ase_results_dir = "/staging/leuven/stg_00096/home/rdewin/ASE/results/") {
  

  # Read the ASE results file for each sample
  ase_resultsfile <- file.path(ase_results_dir, SAMPLEID, paste0(SAMPLEID, "_asereadcounts_pvals_annotated.tsv"))
  ase_results <- read.delim(ase_resultsfile, header = TRUE, as.is = TRUE)
  if (any(grepl(pattern = "chr", x = ase_results$contig))) {
    ase_results$contig <- sub(pattern = "chr", replacement = "", x = ase_results$contig)
  }
  
  # Make results into GRanges object, identify all exonic SNPs and create new df with all of these (contains duplicate SNPs)
  asegr <- GRanges(seqnames = ase_results$contig, ranges = IRanges(start = ase_results$position, end = ase_results$position))
  annothits <- findOverlaps(query = asegr, subject = hsexondb)
  
  # In one case, there were two genes using the same exon ... this just takes the first
  hitgenes <- sapply(mcols(hsexondb[subjectHits(annothits)])$gene_id, FUN = function(x) x[[1]])
  ase_results_annot <- data.frame(ase_results[queryHits(annothits), colnames(ase_results) != "gene"], gene = hitgenes, stringsAsFactors = FALSE)
  
  return(ase_results_annot)
}

combine_pvals_and_adjust <- function(ase_results_annot) {
  outdf <- do.call(rbind, by(data = ase_results_annot, INDICES = ase_results_annot$gene, FUN = combine_pvals))
  outdf$padj <- 1
  outdf[outdf$power, "padj"] <- p.adjust(p = outdf[outdf$power, "pcombined"], method = "fdr")
  return(outdf)
}

add_log2fc_and_gene_names <- function(outdf, l2fcdf, SAMPLEID) {
  outdf[, c("log2fc", "mean_expression", "gene_name")] <- l2fcdf[outdf$gene, c(SAMPLEID, "mean_expression", "gene_name")]
  return(outdf)
}

format_output_data <- function(outdf) {
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  return(outdf)
}

save_output_data <- function(outdf, results_dir, SAMPLEID) {
  # Create the directory if it doesn't exist
  sample_dir <- file.path(results_dir, SAMPLEID)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  
  outfile <- file.path(sample_dir, paste0(SAMPLEID, "_imbalance_expression_vst.txt"))
  write.table(x = outdf, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

save_plot_imbalance <- function(outdf, results_dir, SAMPLEID, significance = 0.01) {
  p1 <- plot_imbalance_expression(imbalancedf = outdf, significance = significance)
  plotfile <- file.path(results_dir, SAMPLEID, paste0(SAMPLEID, "_imbalance_expression_vst_p0.01.png"))
  ggsave(plotfile, p1)
}






