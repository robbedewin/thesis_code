.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")


# Load required libraries
library(ggplot2)
library(GenomicRanges)
library(RColorBrewer)
library(grid)
library(ggrepel)
library(tidyverse)

# Define working directory for ASCAT analysis
ascat_analysis_dir <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/ascat"

# Define ASCAT directory where the segment files are stored
ascat_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat"

# Define output directory
outdir <- "/staging/leuven/stg_00096/home/rdewin/PLOTS/ascat"

# Create output directory if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Load the QC data for plotting
qc_output_file <- file.path(ascat_analysis_dir, "qc_output_file.tsv")
qcdata <- read.delim(qc_output_file, as.is = TRUE)

# Creating a Seqinfo object for CHM13-T2T from genome.dict
genome_dict <- readLines("/staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.dict")

# Parse the genome dictionary data
seq_info_list <- lapply(genome_dict, function(line) {
  if (grepl("^@SQ", line)) {
    parts <- strsplit(line, "\t")[[1]]
    sn <- sub("^SN:(.*)$", "\\1", parts[grep("^SN:", parts)])
    ln <- as.integer(sub("^LN:(.*)$", "\\1", parts[grep("^LN:", parts)]))
    return(list(name = sn, length = ln))
  }
})
# Remove NULL values
seq_info_list <- Filter(Negate(is.null), seq_info_list)

# Create Seqinfo object
seqnames <- sapply(seq_info_list, function(x) x$name)
seqlengths <- sapply(seq_info_list, function(x) x$length)
seqinfoCHM13 <- Seqinfo(seqnames = seqnames, seqlengths = seqlengths, isCircular = rep(FALSE, length(seqnames)))
genome(seqinfoCHM13) <- "CHM13-T2T"

# Read in seg_files and create a GRanges object with annotations
seg_files <- list.files(
  path = ascat_dir,
  pattern = "tumor\\.segments\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

segdata <- setNames(
  lapply(seg_files, function(x) {
    y <- read.delim(x, as.is = TRUE)
    y$chr <- gsub("^([0-9XY]+)$", "chr\\1", as.character(y$chr))
    y$sample <- sub("_tumor", "", y$sample)
    
    gr <- GRanges(
      seqnames = y$chr,
      ranges = IRanges(start = y$startpos, end = y$endpos),
      sample = y$sample,
      nMajor = y$nMajor,
      nMinor = y$nMinor,
      seqinfo = seqinfoCHM13
    )
    gr <- keepSeqlevels(gr, value = paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")
    return(gr)
  }),
  sapply(seg_files, function(x) {
    sub(".*/(P\\d+).*\\.segments\\.txt$", "\\1", x)
  })
)

# Generate a joint genomic region from sample segment data
sampleseg <- unlist(GRangesList(segdata[qcdata$Sample]))
joints <- disjoin(sampleseg)

out <- mapply(sampleid = qcdata$Sample, wgd = qcdata$QC.WGD, nsamples = 1, MoreArgs = list(bins = joints),
              FUN = function(sampleid, wgd, nsamples, bins){
                if (wgd == 0) {
                  gains <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMajor > 1])/nsamples
                  losses <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMinor < 1])/nsamples
                } else {
                  gains <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMajor > 2])/nsamples
                  losses <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMinor < 2])/nsamples
                }
                return(list(gains=gains, losses =losses))
              }, SIMPLIFY = F)

joints$ngains <- Reduce(f = '+', x = lapply(X = out, FUN = function(x) x$gains))
joints$nlosses <- Reduce(f = '+', x = lapply(X = out, FUN = function(x) x$losses))
joints$ncases <- length(qcdata$Sample)
joints$gains <- joints$ngains / joints$ncases
joints$losses <- joints$nlosses / joints$ncases

# Write the joint genomic region to a file
write.table(joints, file = paste0(ascat_analysis_dir, "/ascat_overview.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate cumulative positions
cumdist <- setNames(c(0, cumsum(as.numeric(seqlengths(seqinfoCHM13)[1:22]))), paste0("chr", 1:22))
joints$cumstart <- start(joints) + cumdist[as.character(seqnames(joints))]
joints$cumend <- end(joints) + cumdist[as.character(seqnames(joints))]

# Use names of cumdist as breaks for the plot
breaks_positions <- cumdist[names(cumdist) %in% paste0("chr", 1:22)]

# Create positions for labels on the x-axis in middle of chromosomes
label_position <- breaks_positions[-length(breaks_positions)] + diff(breaks_positions)/2

# Add the last chromosome position manually
last_value <- cumdist[length(cumdist)]
last_chr_midpoint <- last_value + (max(breaks_positions) - last_value)/2
label_position <- c(label_position, last_chr_midpoint)

# Ensure labels_chr has the same length as label_position
labels_chr <- c(1:22)

# Generate ASCAT overview plot
p3 <- ggplot(as.data.frame(joints)) + 
  geom_rect(aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = gains), fill = "red", alpha = .8) +
  geom_rect(aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = -losses), fill = "green", alpha = .8) +
  geom_vline(xintercept = cumdist) + ylim(c(-1, 1)) +
  theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks = label_position, labels = labels_chr) +
  theme(plot.background = element_rect(fill = "white")) +
  labs(x = "", y = "Frequency")

# Save the plot to a file
ggsave(filename = paste0(outdir, "/ascat_overview.png"), plot = p3, width = 10, height = 5)

# Add  labels to position of genes known in T-all to be involved in copy number alterations 
# Load the gene positions
gene_positions_file <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/ascat/genes_position.tsv"
gene_positions <- read.delim(gene_positions_file, header = TRUE, stringsAsFactors = FALSE)

gene_positions$chrom <- paste0("chr", gene_positions$chrom)

# add cumdist to gene_positions
gene_positions$cumdist <- cumdist[gene_positions$chrom] + gene_positions$start


# Add gene positions to the plot
p4 <- p3 + geom_text(data = gene_positions, aes(x = cumdist, y = 0.5, label = gene_id), size = 3, angle = 90, hjust = 0, vjust = 0)

# Save the plot with gene positions to a file
ggsave(filename = paste0(outdir, "/ascat_overview_genes.png"), plot = p4, width = 10, height = 5)


# Convert the GRanges object to a data.frame
joints_df <- as.data.frame(joints)

# Calculate the most aberrated chromosomes
most_aberrated_chromosomes <- joints_df %>%
  group_by(seqnames) %>%
  summarise(gains = sum(gains, na.rm = TRUE), losses = sum(losses, na.rm = TRUE)) %>%
  mutate(aberration_score = gains + losses) %>%
  arrange(desc(aberration_score))

# Print the results
print(most_aberrated_chromosomes)

# Convert the GRanges object to a data.frame
joints_df <- as.data.frame(joints)

# Calculate the length of each chromosome
chrom_lengths <- seqlengths(seqinfoCHM13)

# Add chromosome lengths to the joints dataframe
joints_df <- joints_df %>%
  mutate(chrom_length = chrom_lengths[as.character(seqnames)])

# Calculate the most aberrated chromosomes relative to their length
most_aberrated_chromosomes <- joints_df %>%
  group_by(seqnames) %>%
  summarise(
    gains = sum(gains, na.rm = TRUE),
    losses = sum(losses, na.rm = TRUE),
    chrom_length = first(chrom_length)
  ) %>%
  mutate(
    gains_per_length = gains / chrom_length,
    losses_per_length = losses / chrom_length,
    aberration_score = gains_per_length + losses_per_length
  ) %>%
  arrange(desc(aberration_score))

# Print the results
print(most_aberrated_chromosomes)



# Load cytoband data
cytoband_file <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/ascat/cytoBandMapped.bed"
cytoband <- read.delim(cytoband_file, header = TRUE, col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"))

# Create GRanges object for cytobands
cytoband_gr <- GRanges(
  seqnames = cytoband$chrom,
  ranges = IRanges(start = cytoband$chromStart, end = cytoband$chromEnd),
  name = cytoband$name,
  gieStain = cytoband$gieStain
)

# Define chromosome arms
cytoband_gr$arm <- ifelse(grepl("p", cytoband_gr$name), "p", "q")


# Add cytoband data name to the joints data and also the arm
overlaps <- findOverlaps(joints, cytoband_gr)
joints$cytoband <- NA
joints$arm <- NA
joints$cytoband[subjectHits(overlaps)] <- cytoband_gr$name[queryHits(overlaps)]
joints$arm[subjectHits(overlaps)] <- cytoband_gr$arm[queryHits(overlaps)]




# Function to calculate gains and losses per chromosome arm
calculate_arm_gains_losses <- function(segdata, cytoband_gr) {
  arm_gains_losses <- lapply(names(segdata), function(sample) {
    sample_seg <- segdata[[sample]]
    overlaps <- findOverlaps(sample_seg, cytoband_gr)
    sample_seg$arm <- cytoband_gr$arm[subjectHits(overlaps)]
    
    # Ensure each segment is assigned to a single arm
    sample_seg <- sample_seg[!duplicated(queryHits(overlaps))]
    
    gains <- tapply(sample_seg$nMajor > 1, sample_seg$arm, sum, na.rm = TRUE)
    losses <- tapply(sample_seg$nMinor < 1, sample_seg$arm, sum, na.rm = TRUE)
    
    data.frame(
      sample = sample,
      arm = names(gains),
      gains = gains,
      losses = losses
    )
  })
  
  do.call(rbind, arm_gains_losses)
}

# Calculate gains and losses per chromosome arm
arm_gains_losses <- calculate_arm_gains_losses(segdata, cytoband_gr)


# Save the results to a file
write.table(arm_gains_losses, file.path(outdir, "arm_gains_losses.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Print the results
print(arm_gains_losses)

