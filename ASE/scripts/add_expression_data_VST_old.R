.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

# Load the required libraries
library(DESeq2)
library(GenomicFeatures)
library(ggplot2)

#############################################
# Make counts matrix for DESeq2, normalize and perform variance-stabilizing transformation
# Save normalized counts and log2 fold-change data 
#############################################

# combined counts file from STAR output
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


#############################################
# Combine p-values (of powered SNP loci) per gene
# 
#############################################
## Combine p-values (of powered SNP loci) per gene

source(file = "/staging/leuven/stg_00096/home/rdewin/ASE/scripts/add_expression_data_functions.R")



# Load gene annotations
gtffile <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf"
hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb))) # converts chr1 to 1
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))

# Add in the log2-fold change data and actual gene names
l2fcfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAlog2fc_vst.txt"
l2fcdf <- read.delim(file = l2fcfile, as.is = T)

# Initialize a data frame to store the counts
counts_summary <- data.frame(sample_id = character(), initial_count = integer(), exonic_count = integer(), filtered_count = integer(), stringsAsFactors = FALSE)

# For loop to loop over the sampleIDs
for (SAMPLEID in matchedSamples) {
  #SAMPLEID <- "P011"
  print(paste("Processing sample:", SAMPLEID))
  
  # Read the ASE results file for each sample
  ase_resultsfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_asereadcounts_pvals_annotated.tsv")
  ase_results <- read.delim(ase_resultsfile, header = TRUE, as.is = TRUE)
  if (any(grepl(pattern = "chr", x = ase_results$contig))) {
    ase_results$contig <- sub(pattern = "chr", replacement = "", x = ase_results$contig)
  }

  # Initial count of loci
  initial_count <- nrow(ase_results)

  # make results into GRanges object, identify all exonic SNPs and create new df with all of these (contains duplicate SNPs)
  asegr <- GRanges(seqnames = ase_results$contig, ranges = IRanges(start = ase_results$position, end = ase_results$position))
  annothits <- findOverlaps(query = asegr, subject = hsexondb)
  # converts SNP data into a genomic ranges object, identifies which of these SNPs are located within exonic regions, and stores the overlap information for further analysis. 
  
  # Count after exonic region filtering
  exonic_count <- length(unique(queryHits(annothits)))
  
  # in one case, there were two genes using the same exon ... this just takes the first
  hitgenes <- sapply(mcols(hsexondb[subjectHits(annothits)])$gene_id, FUN = function(x) x[[1]])
  ase_results_annot <- data.frame(ase_results[queryHits(annothits), colnames(ase_results) != "gene"], gene = hitgenes, stringsAsFactors = F)
  
  # create output dataframe with combined p-value per gene + adjust for multiple testing
  outdf <- do.call(rbind, by(data = ase_results_annot, INDICES = ase_results_annot$gene, FUN = combine_pvals))
  
  # Count after filtering for duplicates and power
  filtered_count <- nrow(outdf)
  
  # Add counts to the summary data frame
  counts_summary <- rbind(counts_summary, data.frame(sample_id = SAMPLEID, initial_count = initial_count, exonic_count = exonic_count, filtered_count = filtered_count, stringsAsFactors = FALSE))
  
  outdf$padj <- 1
  outdf[outdf$power, "padj"] <- p.adjust(p = outdf[outdf$power, "pcombined"], method = "fdr")
  
  # add in the log2-fold change data and actual gene names
  outdf[, c("log2fc", "mean_expression", "gene_name")] <- l2fcdf[outdf$gene, c(SAMPLEID, "mean_expression", "gene_name")]

  # format
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  outfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst_new.txt")
  write.table(x = outdf, file = outfile, quote = F, sep = "\t", row.names = F, col.names = T)

  # plot
  p1 <- plot_imbalance_expression(imbalancedf = outdf)
  plotfile <- paste0("/staging/leuven/stg_00096/home/rdewin/ASE/results/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.png")
  ggsave(plotfile, p1)
}

# Write the counts summary to a TSV file
summary_file <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/counts_summary_expression_data.tsv"
write.table(counts_summary, file = summary_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
