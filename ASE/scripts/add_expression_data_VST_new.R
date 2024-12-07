.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

# Load the required libraries
library(DESeq2)
library(GenomicFeatures)
library(ggplot2)

### My combined counts file from STAR output
counts_file <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/counts/combined_counts.tsv"

# Load counts matrix
counts_data <- read.delim(counts_file, header = TRUE, row.names = 1, as.is = TRUE)

# List of matched samples
matchedSamples <- c("P011", "P013", "P016", "P017", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P028", "P029", "P033", "P037", "P041", "P057", "P058", "P059", "P060", "P061", "P064", "P065", "P066", "P086", "P103", "P105")

# Filter counts dataframe
filtered_counts_data <- counts_data[, colnames(counts_data) %in% matchedSamples]

sampleTable <- data.frame(sample = colnames(filtered_counts_data),
                          condition = "T-ALL")

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

# Plot mean expression vs log2 fold-change (to start httpgd server, run `hgd()`)
plot(dds_vstmeans, vst_fc[, 1], main = "Log2 Fold Change vs Mean Expression")
plot(log10(rowMeans(counts(dds, normalized=TRUE))+1), vst_fc[,1])

# Save normalized and VST-transformed data
write.table(assay(dds_vst), 
            file = "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAcounts_vst.txt", 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


# Save log2 fold-change data
l2fcdf <- as.data.frame(vst_fc)

# Create normalized counts dataframe
resdf <- as.data.frame(counts(dds, normalized=TRUE))

# Add gene names (genes were already annotated in the counts file) and mean expression
l2fcdf$gene_name <- rownames(l2fcdf)
l2fcdf$mean_expression <- 2^rowMeans(log2(resdf+1))

# Save log2 fold-change data
write.table(l2fcdf, 
            file = "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAlog2fc_vst.txt", 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Save normalized counts
write.table(resdf, 
            file = "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAcounts_normalised_T-ALL.txt", 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)



## Combine p-values (of powered SNP loci) per gene

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

plot_imbalance_expression <- function(imbalancedf) {
  imbalancedf <-  imbalancedf[order(imbalancedf$mean_expression, decreasing = F), ]
  labeldf <- data.frame(pos = unlist(lapply(X = 10^(0:4), FUN = function(x) sum(imbalancedf$mean_expression < x))), expr = 10^(0:4), stringsAsFactors = F)
  
  outdf_bak <- imbalancedf
  imbalancedf <- imbalancedf[!grepl(pattern = "^HLA.*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^IG[HLK].*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^TR[ABDG][VCDJ].*", x = imbalancedf$gene_name, perl = T), ]

  imbalancedf$notes <- ifelse(imbalancedf$padj > 0.05, "nonsig", 
                              ifelse(imbalancedf$log2fc >= 1, "up",
                                     ifelse(imbalancedf$log2fc <= -.73, "down", "nonsig")))
  
  p1 <- ggplot(data = imbalancedf, mapping = aes(x = 1:nrow(imbalancedf), y = -sign(log2fc)*log10(pcombined)))
  p1 <- p1 + geom_point(mapping = aes(colour = notes, size = abs(log2fc)), 
                        show.legend = F, alpha = .4)
  p1 <- p1 + geom_hline(yintercept = c(-1,1)*-log10(max(imbalancedf[imbalancedf$padj < .05, "pcombined"])), linetype = "dashed", colour = "grey") +
    geom_text(data = imbalancedf[imbalancedf$notes != "nonsig", ], mapping = aes(x = which(imbalancedf$notes != "nonsig"), y = -sign(log2fc)*log10(pcombined), label = gene_name), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(imbalancedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_y_continuous(breaks = seq(-10,10,2), oob = scales::squish, limits = c(-10,10))
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$expr, name = "mean expression (normalised)")
  # p1 <- p1 + scale_color_brewer(type = "div", palette = "RdBu", direction = -1)
  p1 <- p1 + scale_color_manual(values = c(nonsig = "#e0e0e0", up = "#ef8a62", down = "#67a9cf"))
  p1 <- p1 + scale_size_continuous(range = c(1,7.5))
  p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = -90)) + labs(x = NULL)
  return(p1)  
}


# Load gene annotations
gtffile <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf"
hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb))) # converts chr1 to 1
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))

# Add in the log2-fold change data and actual gene names
l2fcfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAlog2fc_vst.txt"
l2fcdf <- read.delim(file = l2fcfile, as.is = T)

# For loop to loop over the sampleIDs
for (SAMPLEID in matchedSamples) {
  #SAMPLEID <- "P013"
  print(paste("Processing sample:", SAMPLEID))
  
  # Read the ASE results file for each sample
  ase_resultsfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_asereadcounts_nomatch_pvals_annotated.tsv")
  ase_results <- read.delim(ase_resultsfile, header = TRUE, as.is = TRUE)
  if (any(grepl(pattern = "chr", x = ase_results$contig))) {
    ase_results$contig <- sub(pattern = "chr", replacement = "", x = ase_results$contig)
  }
  
  # make results into GRanges object, identify all exonic SNPs and create new df with all of these (contains duplicate SNPs)
  asegr <- GRanges(seqnames = ase_results$contig, ranges = IRanges(start = ase_results$position, end = ase_results$position))
  annothits <- findOverlaps(query = asegr, subject = hsexondb)
  # converts SNP data into a genomic ranges object, identifies which of these SNPs are located within exonic regions, and stores the overlap information for further analysis. 
  
  # in one case, there were two genes using the same exon ... this just takes the first
  hitgenes <- sapply(mcols(hsexondb[subjectHits(annothits)])$gene_id, FUN = function(x) x[[1]])
  ase_results_annot <- data.frame(ase_results[queryHits(annothits), colnames(ase_results) != "gene"], gene = hitgenes, stringsAsFactors = F)
  
  # create output dataframe with combined p-value per gene + adjust for multiple testing
  outdf <- do.call(rbind, by(data = ase_results_annot, INDICES = ase_results_annot$gene, FUN = combine_pvals))
  outdf$padj <- 1
  outdf[outdf$power, "padj"] <- p.adjust(p = outdf[outdf$power, "pcombined"], method = "fdr")
  
  # add in the log2-fold change data and actual gene names
  outdf[, c("log2fc", "mean_expression", "gene_name")] <- l2fcdf[outdf$gene, c(SAMPLEID, "mean_expression", "gene_name")]

  # format
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  outfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.txt")
  write.table(x = outdf, file = outfile, quote = F, sep = "\t", row.names = F, col.names = T)

  # plot
  p1 <- plot_imbalance_expression(imbalancedf = outdf)
  plotfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.png")
  ggsave(plotfile, p1)
    
}

