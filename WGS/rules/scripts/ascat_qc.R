.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")
#library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(GenomicRanges)
library(RColorBrewer)
library(parallel)
library(grid)
library(ggrepel)
library(stringr)
library(tidyverse)

# Open the log file
#log <- file(snakemake@log[[1]], open = "wt")
#sink(log)
#sink(log, type = "message")

# Define input files and output directories based on snakemake parameters
#sample_ids <- snakemake@params[["sample_ids"]]
#ascatfiles <- dirname(snakemake@params[["ascat_dir"]])

#outdir <- dirname(snakemake@output[[1]])
#qc_output_file <- snakemake@output[[1]]
#seg_plot <- snakemake@output[[2]]

# Hard-coded sample IDs and paths for debugging or standalone runs
sample_ids <- c("P028", "P022", "P024", "P018", "P019", "P020", "P023", "P026", "P013", "P016")
ascatfiles <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/"

outdir <- dirname("/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/qc/summary.csv")
qc_output_file <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/qc/summary.csv"
seg_plot <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/qc/segments_plot.png"

# Create output directory if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

setwd(ascatfiles)
getwd()

# Store all segmentsfiles in a list
seg_files <- list.files(
  path = file.path(ascatfiles, list.dirs(path = ascatfiles, full.names = FALSE)[list.dirs(path = ascatfiles, full.names = FALSE) %in% sample_ids]),
  pattern = ".segments.txt$",
  full.names = TRUE,
  recursive = TRUE
)
seg_files
# Store all qc files in a list
qc_files <- list.files(
  path = file.path(ascatfiles, list.dirs(path = ascatfiles, full.names = FALSE)[list.dirs(path = ascatfiles, full.names = FALSE) %in% sample_ids]),
  pattern = "_objects.Rdata$",
  full.names = TRUE,
  recursive = TRUE
)
qc_files


qc_data <- lapply(qc_files, function(x) {
  load(x)
  data.frame(Sample = sub(".*/(P\\d+)_ASCAT_objects\\.Rdata$", "\\1", x), QC = QC)
}) %>% bind_rows()
qc_data

# Write the qc_data data frame to a tab-separated .csv file
write.table(qc_data, qc_output_file, sep = "\t", row.names = FALSE, quote = FALSE)

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
print(seqinfoCHM13)


# Read in seg_files and create a GRanges object with annotations
segdata <- setNames(  # Set names for each GRanges object in the list
  lapply(seg_files, function(x) {
    y <- read.delim(x, as.is = TRUE)
    # Ensure chromosome names are correctly prefixed with 'chr' if necessary
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
  # Extract the sample IDs from the filenames to use as names for GRanges objects
  sapply(seg_files, function(x) {
    sub(".*/(P\\d+).*\\.segments\\.txt$", "\\1", x)
  })
)


# Load the QC data for plotting
qcdata <- read.delim(qc_output_file, as.is = TRUE)

#Plotting loh and ploidy
m <- (2.9 - 1)/(0-0.93) #??

# Define and create the plot
p1 <- ggplot(data = qcdata, mapping = aes(x = QC.LOH, y = QC.ploidy, colour = QC.WGD))
p1 <- p1 + geom_point() +
  geom_abline(slope = m, intercept = 2.9) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

# Add labels with geom_label_repel
p1 <- p1 + geom_label_repel(
  aes(label = Sample),  
  label.size = 0.3,     # Controls the border size around labels, set to NA to remove
  size = 1.8,           # Text size
  max.overlaps = Inf    # Allows for infinite overlaps; adjust as necessary
)

ggsave(filename = "qc/annotated_frac_loh_vs_ploidy.png", plot = p1, width = 10, height = 7)

# Plot LOH vs Ploidy
plotLOH <-  ggplot(qcdata, aes(x = QC.LOH, y = QC.ploidy, color = QC.WGD)) +
              geom_point() +
              geom_abline(slope = (2.9 - 1)/(0 - 0.93), intercept = 2.9) +
              theme_minimal() +
              theme(plot.background = element_rect(fill = "white")) +
              geom_label_repel(aes(label = Sample), label.size = 0.3, size = 1.8, max.overlaps = Inf) +
              ggsave("qc/fraction_loh_vs_ploidy.png", width = 10, height = 7)

#Fraction of genome altered
# Check if segment files are non-empty and filter out empty files
checkfiles <- sapply(seg_files, FUN = function(x) file.size(x) > 1)
seg_files <- seg_files[checkfiles]

# Calculate the fractional aberration for diploid and tetraploid states, and LOH
fracaberdip <- sapply(segdata, FUN = function(x) sum(width(x[x$nMajor != 1 & x$nMinor != 1])) / sum(width(x)))
fracabertet <- sapply(segdata, FUN = function(x) sum(width(x[x$nMajor != 2 & x$nMinor != 2])) / sum(width(x)))
fracloh <- sapply(segdata, FUN = function(x) sum(width(x[xor(x$nMinor == 0, x$nMajor == 0)])) / sum(width(x)))


# Organize output into a data frame
outdf1 <- data.frame(
  sampleid = sub("_tumor.segments.txt", "", basename(seg_files)),
  frac_aber_dip = fracaberdip,
  frac_aber_tet = fracabertet,
  frac_loh = fracloh
)

print(outdf1)


#Integrating Aberration Data into QC Data
qcdata$QC.frac_aber <- ifelse(qcdata$QC.WGD, outdf1$frac_aber_tet, outdf1$frac_aber_dip)
qcdata$QC.frac_aber[is.na(qcdata$ploidy)] <- NA
qcdata

# Identify samples with specific characteristics, e.g., presence of MYC amplifications
MYC_pos <- sample_ids
qcdata$MYC_sig <- ifelse(qcdata$Sample %in% MYC_pos, 1, 0)

# Prepare a color palette
colpal <- colorRampPalette(brewer.pal(12, "Set2"))
myPal <- colpal(length(unique(qcdata$Sample)))

# Plot fractional aberration with indication of MYC signaling
p1 <- ggplot(qcdata, aes(x = Sample, y = QC.frac_aber, size = 0.2, color = MYC_sig)) + 
  geom_jitter(width = 0.25, show.legend = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("qc/Frac_aberrated_with_Myc.png", plot = p1, width = 10, height = 5)

# Plot genetic instability index
p2 <- ggplot(qcdata, aes(x = Sample, y = QC.GI, size = 0.2, color = MYC_sig)) + 
  geom_jitter(width = 0.25, show.legend = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("qc/QC.GI_of_the_cohort.png", plot = p2, width = 10, height = 5)


p2 <- ggplot(data = qcdata, mapping = aes(x = Sample, y = QC.GI, size = 0.2)) + 
  geom_jitter(mapping = aes(color = Sample), width = 0.25, show.legend = TRUE)
p2 <- p2 + theme(axis.text.x = element_text(size = 12))  # Adjust the size (e.g., size = 12)
p2 <- p2 + scale_fill_manual(values = myPal) + scale_color_manual(values = myPal) + theme_minimal() + theme(plot.background = element_rect(fill = "white")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p2
ggsave(filename = "qc/QC.GI.png", plot = p2, width = 10, height = 5)

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
joints
joints$gains <- joints$ngains / joints$ncases
joints$losses <- joints$nlosses / joints$ncases

# Calculate cumulative positions
cumdist <- setNames(c(0, cumsum(as.numeric(seqlengths(seqinfoCHM13)[1:22]))), paste0("chr", 1:22))
joints$cumstart <- start(joints) + cumdist[as.character(seqnames(joints))]
joints$cumend <- end(joints) + cumdist[as.character(seqnames(joints))]

# Use names of cumdist as breaks for the plot
breaks_positions <- cumdist[names(cumdist) %in% paste0("chr", 1:22)]
labels_chr <- 1:22

# Generate ASCAT overview plot
# Generate the ASCAT overview plot
p3 <- ggplot(as.data.frame(joints)) + 
  geom_rect(aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = gains), fill = "#fc8d59", alpha = .6) +
  geom_rect(aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = -losses), fill = "#91bfdb", alpha = .6) +
  geom_vline(xintercept = cumdist) + ylim(c(-1, 1)) +
  theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks = breaks_positions, labels = labels_chr) +
  theme(plot.background = element_rect(fill = "white")) +
  labs(x = "", y = "Frequency")
ggsave("qc/ascat_overview.png", plot = p3, width = 10, height = 5)


#Additional QC from the ASCAT plots
#QC was determined as follows:
#  - **No ASCAT fit**: ASCAT could not find an optimal ploidy and purity value.
#- **Likely normal**: CNA profile matches with a normal sample. Filter: LOH<0.1 & GI<0.1 & purity=1 & ((sex='XX' & ploidy>=1.99 & ploidy<=2.01) | (sex='XY' & ploidy>=1.945 & ploidy<=1.965)).
#- **Noisy logR**: logR track (either tumour or normal) is noisy. Filter: tumour_mapd>=0.75 | normal_mapd>=0.75.
#- **Likely wrong germline**: matched normal is likely to come from another patient. Filter: tumour_mapd>0.4 & normal_mapd>0.4 & frac_homo>0.1.
#- **Oversegmented**: difference in segmentation between logR and BAF tracks. Filter: n_segs_logRBAF_diff>250.
#- **Wrong fit**: CNA profile is incorrect, with extreme losses. Filter: mode_majA=0.
#- **HD size**: CNA profile contains large segments with homozygous deletion. Filter: homdel_largest>=20e6 | homdel_size>=40e6.
#- **Contamination or swap** and **Large CNV**: germline data contains copy-number changes. Filter: we ran aspcf on all germline samples and manually reviewed cases where events were spotted. Cases with large stretches of homozygosity were considered as valid and were not discarded.
#- **Pass**: sample has correct metrics.

# Applying QC conditions based on ASCAT analysis outputs
qcdata$Likely_normal <- ifelse(
  qcdata$QC.LOH < 0.1 & qcdata$QC.GI < 0.1 & qcdata$QC.purity == 1 & (
    (qcdata$QC.sex == 'XX' & qcdata$QC.ploidy >= 1.99 & qcdata$QC.ploidy <= 2.01) |
    (qcdata$QC.sex == 'XY' & qcdata$QC.ploidy >= 1.945 & qcdata$QC.ploidy <= 1.965)
  ), TRUE, FALSE
)

qcdata$Noisy_LogR <- ifelse(
  qcdata$QC.tumour_mapd >= 0.75 | (!is.na(qcdata$QC.normal_mapd) & qcdata$QC.normal_mapd >= 0.75),
  TRUE, FALSE
)

qcdata$Likely_wrong_Germline <- ifelse(
  qcdata$QC.tumour_mapd > 0.4 | (!is.na(qcdata$QC.normal_mapd) & qcdata$QC.normal_mapd >= 0.4) & qcdata$QC.frac_homo > 0.1,
  TRUE, FALSE
)

qcdata$Oversegmented <- ifelse(
  qcdata$QC.n_segs_logR - qcdata$QC.n_segs_BAF > 250,
  TRUE, FALSE
)

qcdata$wrong_fit <- ifelse(qcdata$QC.mode_majA == 0, TRUE, FALSE)

qcdata$HD_size <- ifelse(
  qcdata$QC.homdel_largest >= 20e6 | qcdata$QC.homdel_size >= 40e6,
  TRUE, FALSE
)

# Write updated QC data to a file
write.table(qcdata, "qc_data.csv", sep = "\t", row.names = FALSE, quote = FALSE)


m <- (2.9 - 1)/(0-0.93)
# Plotting segmentation differences
p1 <- ggplot(qcdata, aes(x = QC.n_segs_BAF, y = QC.n_segs_logR, colour = QC.WGD)) +
  geom_point() +
  geom_abline(slope = m, intercept = 250) + 
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white")) +
  geom_label_repel(aes(label = Sample), label.size = 0.3, size = 1.8, max.overlaps = Inf)

ggsave("qc/oversegmented.pdf", plot = p1, width = 10, height = 5)

# Making Heatmap of the CNV landscape (following some steps of cnSpec)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
library(circlize)

# Define genomic ranges and create windows of 5 million bases
chr_df <- read.chromInfo()$df
chr_df <- chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_window <- makeWindows(chr_gr, w = 5e6)

# Define a new window based on joint data
new_window <- GRanges(seqnames = seqnames(joints), ranges = ranges(joints))

# Define a function to calculate average signals over windows, handling both numeric and character vectors
average_in_window <- function(window, gr, v, empty_v = NA, method = "weighted") {
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    }
  }
  
  return(u)
}
sampleseg$total_CN <- sampleseg$nMajor + sampleseg$nMinor
PTCL1_test <- sampleseg[sampleseg$sample == "PTCL1_T",]
#sampleseg
#CN_data <- average_in_window(new_window,PTCL1_test,PTCL1_test$total_CN)

# Calculate total copy number and apply the function to segment data
CNsegdata <- lapply(segdata, function(x) {
  x$total_CN <- x$nMinor + x$nMajor
  x
})


total_CN_data <- lapply(CNsegdata, function(x) {
  average_in_window(new_window, x, x$total_CN)
})

# Prepare data frame for Complex Heatmap
chr_window_df <- as.data.frame(new_window)
sample_names <- names(total_CN_data)

# Populate the data frame with CN data for each sample
for (sample_name in sample_names) {
  chr_window_df[[sample_name]] <- total_CN_data[[sample_name]]
}



final_table <- chr_window_df %>% select(-width, -strand)
long_table <- final_table %>% pivot_longer(cols = starts_with("P"), names_to = "sample", values_to = "cn")
long_table <- long_table %>% rename(chromosome = "seqnames", segmean = "cn")

# Define a color function for heatmap and generate heatmap
col_fun <- colorRamp2(c(0, 2, 6), c("blue", "white", "red"))
p1 <- cnSpec(long_table, y = genomeBoundaries, genome = "CHM13-T2T", CNscale = "absolute")
ggsave("heatmap_first_time.pdf", plot = p1, width = 10, height = 5)

# Assuming long_table has columns named 'chromosome', 'start', 'end', 'sample', and 'cn' for copy number
# Creating a matrix from the long_table for heatmap plotting
heatmap_matrix <- with(long_table, matrix(cn, nrow = length(unique(chromosome)), ncol = length(unique(sample)), dimnames = list(unique(chromosome), unique(sample))))

# Create the heatmap
heatmap <- Heatmap(heatmap_matrix, name = "CNV", col = colorRamp2(c(0, 2, 6), c("blue", "white", "red")), 
                   row_title = "Chromosome", column_title = "Sample", 
                   cluster_rows = FALSE, cluster_columns = FALSE)

draw(heatmap)


#GenVisR heatmap
BiocManager::install("GenVisR")
library(GenVisR)

genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg38",], max)
genomeBoundaries$chromStart <- 0
colnames(genomeBoundaries) <- c("chromosome", "end", "start")
cumdist

replace_na_with_2 <- function(x) {
  x[is.na(x)] <- 2
  return(x)
}
replace_na_with_2(long_table$segmean)
p1 <- cnSpec(long_table, y = genomeBoundaries, genome = "hg38", CNscale = "absolute")
p1
ggsave(filename = "heatmap_first_time.pdf",plot = p1, width =10, height= 5 )
setwd("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/plots")
genomeBoundaries_2 <- genomeBoundaries
genomeBoundaries_2$start <- genomeBoundaries_2$start - 1e8
genomeBoundaries_2$end <- genomeBoundaries_2$end + 1e8
