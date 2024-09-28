# Introduction
# Allele-Specific Expression (ASE) Analysis Pipeline
# This notebook provides a comprehensive pipeline for performing Allele-Specific Expression (ASE) analysis. 
# The pipeline processes raw sequencing data to identify loci with significant allele-specific expression patterns. 
# It includes filtering allele counts, combining loci information, running ASEReadCounter, computing statistical significance, 
# annotating results, and visualizing findings through a Manhattan plot.

# Goals:
# 
# Filter and prepare allele counts for ASE analysis.
# Combine loci information and generate VCF files.
# Utilize GATK's ASEReadCounter for read count data.
# Compute p-values to identify significant allele imbalance.
# Annotate significant loci with gene information.
# Visualize results using a Manhattan plot with annotations.

# Inputs:
# 
# Reference Alleles Files: Contain allele information for each chromosome.
# Sample Allele Counts Files: Allele counts from sequencing data.
# RNA BAM Files: Aligned RNA sequencing reads.
# Reference Genome Fasta File: The reference genome sequence.
# Gene Annotation GTF File: Contains gene annotations for the genome.

# Outputs:
# 
# Filtered Allele Counts Files: Allele counts after applying filters.
# Combined Loci Information: Text and VCF files with combined allele counts.
# ASE Read Counts: Output from ASEReadCounter.
# Annotated ASE Results: Data frame with computed p-values and gene annotations.
# Manhattan Plot: Visualization of significant loci across the genome.

# Set library path
.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")


# Load necessary libraries
library(readr)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(BSgenome.Hsapiens.UCSC.hs1)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(VGAM)
library(BiocGenerics)
library(parallel)
library(S4Vectors)


# Source the utility functions
source(file = "/staging/leuven/stg_00096/home/rdewin/WGS/rules/scripts/utils.R") 



# Set variables
reference_alleles_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/ASCAT/ReferenceFiles"
ascat_counts_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat"
alias <- "tumor"
min_depth <- 5
results_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/results"
RNA_dir <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/star"
WGS_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources"
gatk_exe <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/gatk"
java_cmd <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/java"
gatk_jar <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar"
java_home <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS"
gtf_file <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf"
filter_cutoff <- 0.01
sig_threshold <- -log10(0.05)

#Get list of samples
rna_files <- list.files(RNA_dir, pattern = "_Aligned.sortedByCoord.out.bam$", full.names = TRUE, recursive = TRUE)
rna_samples <- basename(dirname(rna_files))

dna_files <- list.dirs(ascat_counts_dir, full.names = TRUE, recursive = FALSE)
dna_samples <- basename(dna_files)

common_samples <- intersect(rna_samples, dna_samples)

sample_id <- common_samples[1]

# Ensure output directories exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Define chromosomes to process
chromosomes <- c(1:22, "X")

get_allele_counts_filtered <- function(chr, reference_alleles_dir, ascat_counts_dir, min_depth = 3, sample_id, alias = "tumor", output_dir) {
  library(readr)
  
  # Conversion map from numeric to nucleotide
  allele_conversion <- c("1" = "A", "2" = "C", "3" = "G", "4" = "T")
  
  # Construct file paths
  reference_alleles_file <- file.path(reference_alleles_dir, paste0("allele_T2T_chr", chr, ".txt"))
  sample_allele_counts_file <- file.path(ascat_counts_dir, sample_id, "alleleFrequencies", paste0(sample_id, "_", alias, "_alleleFrequencies_chr", chr, ".txt"))
  
  # Check if input files exist
  if (!file.exists(reference_alleles_file)) {
    warning(paste("Reference alleles file not found:", reference_alleles_file))
    return(NULL)
  }
  if (!file.exists(sample_allele_counts_file)) {
    warning(paste("Sample allele counts file not found:", sample_allele_counts_file))
    return(NULL)
  }
  
  # Output file path
  output_file <- file.path(output_dir, "allele_counts", paste0(sample_id, "_", alias, "_filtered_allele_counts_chr", chr, ".txt"))
  
  # Read reference alleles
  reference_alleles <- read_tsv(
    file = reference_alleles_file,
    col_names = TRUE,
    col_types = "iii"
  )
  
  colnames(reference_alleles) <- c("pos", "ref_num", "alt_num")

  # Convert numeric alleles to nucleotides
  reference_alleles <- reference_alleles %>%
    mutate(
      ref = allele_conversion[as.character(ref_num)],
      alt = allele_conversion[as.character(alt_num)]
    ) %>%
    select(pos, ref, alt)
  
  # Read sample allele counts
  sample_allele_counts <- read_tsv(
    file = sample_allele_counts_file,
    col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"),
    col_types = "ciiiiii",
    comment = "#"
  )
  
  # Combine reference alleles and sample counts
  allele_counts <- left_join(reference_alleles, sample_allele_counts, by = "pos")
  
  

  # Calculate counts for ref and alt alleles
  allele_counts <- allele_counts %>%
    mutate(
      count_ref = case_when(
        ref == "A" ~ count_A,
        ref == "C" ~ count_C,
        ref == "G" ~ count_G,
        ref == "T" ~ count_T
      ),
      count_alt = case_when(
        alt == "A" ~ count_A,
        alt == "C" ~ count_C,
        alt == "G" ~ count_G,
        alt == "T" ~ count_T
      )
    )
  
  allele_counts <- allele_counts %>%
    select(chr, pos, ref, alt, count_ref, count_alt)

  # Filter based on minimum depth
  filtered_counts <- allele_counts %>%
    filter(count_ref >= min_depth, count_alt >= min_depth)
  
  if (nrow(filtered_counts) == 0) {
    warning(paste("No loci passed the minimum depth filter for sample", sample_id, "chr", chr))
    return(NULL)
  }
  
  # Write filtered counts to output file
  write_tsv(filtered_counts, output_file, col_names = FALSE)
  
  message(paste("Filtered allele counts written to:", output_file))
  
  return(NULL)
}

combine_loci_nomatch <- function(sample_id, results_dir) {
  library(readr)
  library(VariantAnnotation)
  library(GenomicRanges)
  
  # Define input and output file paths
  allele_counts_files <- list.files(
    path = file.path(results_dir, sample_id, "allele_counts"),
    pattern = paste0(alias, "_filtered_allele_counts_chr.*\\.txt$"),
    full.names = TRUE
  )
  
  combined_counts_file <- file.path(results_dir, sample_id, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  combined_vcf_file <- file.path(results_dir, sample_id, paste0(sample_id, "_hetSNPs_nomatch.vcf"))
  
  # Read and combine allele count files
  combined_counts <- tryCatch(
    {
      combined_data <- do.call(
        rbind,
        lapply(
          X = allele_counts_files,
          FUN = function(file) {
            read_tsv(
              file,
              col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"),
              col_types = "ciccii",
              progress = FALSE
            )
          }
        )
      )
      if (nrow(combined_data) == 0) {
        stop("Combined allele count data is empty.")
      }
      combined_data
    },
    error = function(e) {
      stop(paste("Error reading allele count files:", e$message))
    }
  )


  # Write combined counts to file
  write_tsv(combined_counts, combined_counts_file, col_names = TRUE)
  
  # Create VCF object
  gr <- GRanges(seqnames = combined_counts$chr, ranges = IRanges(start = combined_counts$pos, width = 1))
  
  # Set seqlevelsStyle and seqinfo
  seqlevelsStyle(gr) <- seqlevelsStyle(BSgenome.Hsapiens.UCSC.hs1)
  seqinfo(gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hs1)[seqlevels(gr)]

  # Create REF and ALT columns
  ref <- DNAStringSet(combined_counts$ref)
  alt <- CharacterList(as.list(combined_counts$alt)) 

  # Create QUAL and FILTER
  qual <- rep(NA_real_, nrow(combined_counts))
  filter <- rep("PASS", nrow(combined_counts))
  
  # Create fixed DataFrame
  fixed <- DataFrame(REF = ref, ALT = alt, QUAL = qual, FILTER = filter)
  
  # Create genotype matrices
  gt <- matrix("0/1", nrow = nrow(combined_counts), ncol = 1)
  colnames(gt) <- sample_id
  
  ad_array <- array(dim = c(nrow(combined_counts), 1, 2))
  ad_array[,,1] <- as.integer(combined_counts$count_ref)
  ad_array[,,2] <- as.integer(combined_counts$count_alt)
  dimnames(ad_array) <- list(NULL, sample_id, c("Ref", "Alt"))
  
  dp <- matrix(as.integer(combined_counts$count_ref + combined_counts$count_alt), nrow = nrow(combined_counts), ncol = 1)
  colnames(dp) <- sample_id
  
  # Create geno SimpleList
  geno_list <- SimpleList(GT = gt, AD = ad_array, DP = dp)
  
  # Create colData
  col_data <- DataFrame(row.names = sample_id)
  
  # Create VCF object
  vcf_obj <- VCF(
    rowRanges = gr,
    colData = col_data,
    fixed = fixed,
    geno = geno_list
  )
  
  # Create meta-information lines without 'fileDate' and 'fileformat'
  meta_info <- DataFrame(
    row.names = c("phasing", "source"),
    Value = c(
      "unphased",
      paste("VariantAnnotation", packageVersion("VariantAnnotation"))
    )
  )
  
  # Create FORMAT header
  format_df <- DataFrame(
    Number = c("1", "2", "1"),
    Type = c("String", "Integer", "Integer"),
    Description = c(
      "Genotype",
      "Allelic depths (number of reads in each observed allele)",
      "Total read depth"
    ),
    row.names = c("GT", "AD", "DP")
  )
  
  # Create VCFHeader
  vcf_header <- VCFHeader(
    samples = sample_id,
    header = DataFrameList(
      META = meta_info,
      FORMAT = format_df
    )
  )
  
  # Set header in VCF object
  header(vcf_obj) <- vcf_header
  
  # Write VCF file
  writeVcf(vcf_obj, combined_vcf_file , index = FALSE)
  
  # Apply sed transformation: delete existing ##fileformat and add new line at the top
  system(paste("sed '/^##fileformat/ d' ", combined_vcf_file, " | sed '1i ##fileformat=VCFv4.3' > tmp_vcf && mv tmp_vcf ", combined_vcf_file))

  # zip and index the VCF file using bgzip and tabix (need to install these tools in conda)
  system(paste("bgzip -f ", combined_vcf_file))
  system(paste("tabix -p vcf ", paste0(combined_vcf_file, ".gz")))
  
  message(paste("Combined VCF written to:", paste0(combined_vcf_file, ".gz")))
}


run_ASEReadCounter <- function(sample_id, results_dir, RNA_dir, WGS_dir, gatk_jar, java_cmd, java_home, min_base_quality = 20, min_mapping_quality = 35) {
  # Set JAVA_HOME environment variable
  Sys.setenv(JAVA_HOME = java_home)
  
  # Define file paths
  het_snp_vcf <- file.path(results_dir, sample_id, paste0(sample_id, "_hetSNPs_nomatch.vcf.gz"))
  bam_file <- file.path(RNA_dir, sample_id, paste0(sample_id, "_Aligned.sortedByCoord.withRG.bam"))
  ref_genome <- file.path(WGS_dir, "genome.fa")
  
  output_dir <- file.path(results_dir, sample_id)
  output_file <- file.path(output_dir, paste0(sample_id, "_asereadcounts_nomatch.tsv"))
  
  # Construct the command
  cmd <- paste(
    shQuote(java_cmd), "-Xmx32G -jar", shQuote(gatk_jar), "ASEReadCounter",
    "-I", shQuote(bam_file),
    "-V", shQuote(het_snp_vcf),
    "-R", shQuote(ref_genome),
    "-O", shQuote(output_file),
    "--min-depth 1",
    "--min-mapping-quality", min_mapping_quality,
    "--min-base-quality", min_base_quality,
    "--lenient",
    " 2> ", shQuote(file.path(output_dir, "ASEReadCounter.stderr"))
  )
  
  # Record and display the start time
  start_time <- Sys.time()
  message("Starting ASEReadCounter execution.", start_time)

  # Execute the command
  system(cmd)

  # Record the end time
  end_time <- Sys.time()
  message("ASEReadCounter execution completed.", end_time)
     
  message(paste("ASEReadCounter output written to:", output_file))
}


compute_pvals <- function(sample_id, results_dir, filter_cutoff = 0.01) {
  library(readr)
  library(VGAM)
  
  # Define file paths
  ase_counts_file <- file.path(results_dir, sample_id, paste0(sample_id, "_asereadcounts_nomatch.tsv"))
  genome_counts_file <- file.path(results_dir, sample_id, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  
  # Read counts data
  ase_counts <- read_tsv(file = ase_counts_file, col_names = TRUE, col_types = "ciccciiiiiiii")
  
  genome_counts <- read_tsv(file = genome_counts_file, col_names = TRUE, col_types = "ciccii")

  colnames(genome_counts) <- c("chr", "pos", "ref", "alt", "refCountGenome", "altCountGenome")

  # Merge data
  asedf <- left_join(ase_counts, genome_counts, by = c("contig" = "chr", "position" = "pos"))
  
  # Select relevant columns
  asedf <- asedf %>%
    select(contig, position, refAllele, altAllele, refCountGenome, altCountGenome, refCount, altCount)

  # Compute filter and p-values
  asedf <- asedf %>%
    rowwise() %>%
    mutate(
      filter = qbeta(c(filter_cutoff / 2, 1 - filter_cutoff / 2), shape1 = refCountGenome + 1, shape2 = altCountGenome + 1)[1],
      pval = betabinom.test.ab(
        q = refCount, size = refCount + altCount,
        shape1 = refCountGenome + 1, shape2 = altCountGenome + 1, alternative = "two.sided"
      ),
      padj = p.adjust(pval, method = "fdr")
    )
  
  # Write results to file
  output_file <- file.path(results_dir, sample_id, paste0(sample_id, "_asereadcounts_nomatch_pvals.tsv"))
  write_tsv(asedf, output_file, col_names = TRUE)

  return(asedf)
}

annotate_ase_results <- function(asedf, gtf_file) {
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  
  # Read gene annotations
  gtf <- rtracklayer::import(gtf_file)
  transcripts_gtf <- gtf[gtf$type == "transcript"]
  
  # Create GRanges for ASE loci
  ase_gr <- GRanges(
    seqnames = asedf$contig,
    ranges = IRanges(start = asedf$position, end = asedf$position)
  )
  
  # Harmonize seqlevels
  seqlevelsStyle(ase_gr) <- seqlevelsStyle(transcripts_gtf)
  
  # Find overlaps
  overlaps <- findOverlaps(query = ase_gr, subject = transcripts_gtf)
  
  # Convert overlaps to data.table
  overlap_dt <- data.table(
    query = queryHits(overlaps),
    subject = subjectHits(overlaps)
  )
  
  # Add gene names
  overlap_dt[, gene := mcols(transcripts_gtf)$gene_name[subject]]
  
  # Collapse gene names per query
  gene_list <- overlap_dt[, .(gene = paste(unique(gene), collapse = ",")), by = query]
  
  # Assign gene names to asedf
  asedf$gene <- NA
  asedf$gene[gene_list$query] <- gene_list$gene

  # Write annotated results to file
  output_file <- file.path(results_dir, sample_id, paste0(sample_id, "_asereadcounts_nomatch_pvals_annotated.tsv"))
  write_tsv(asedf, output_file, col_names = TRUE)
  
  return(asedf)
}

plot_ase_manhattan <- function(asedf, sig_threshold = -log10(0.05)) {
  library(ggplot2)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hs1)
  
  # Remove 'chr' prefix
  asedf$contig <- gsub("^chr", "", asedf$contig)

  # Ensure contig is a factor with ordered levels
  chroms <- c(as.character(1:22), "X")
  asedf$contig <- factor(asedf$contig, levels = chroms)

  # Create a mapping of chromosomes to numeric indices
  chrom_map <- data.frame(contig = chroms, chrom_index = 1:length(chroms), stringsAsFactors = FALSE)

  # Merge chrom_index into asedf
  asedf <- left_join(asedf, chrom_map, by = "contig")
  
  # Get chromosome lengths from BSgenome
  bsgenome_hs1 <- BSgenome.Hsapiens.UCSC.hs1
  chrom_lengths <- seqlengths(bsgenome_hs1)
  
  # Keep only chromosomes of interest
  chroms_full <- paste0("chr", chroms)
  chrom_lengths <- chrom_lengths[names(chrom_lengths) %in% chroms_full]
  
  # Remove 'chr' prefix from names
  names(chrom_lengths) <- gsub("^chr", "", names(chrom_lengths))
  
  # Create data frame
  chrom_lengths_df <- data.frame(
    contig = names(chrom_lengths),
    chr_len = as.numeric(chrom_lengths),
    stringsAsFactors = FALSE
  )
  
  # Add chrom_index
  chrom_lengths_df <- left_join(chrom_lengths_df, chrom_map, by = "contig")
  
  # Arrange and compute chromosome offsets
  chrom_lengths_df <- chrom_lengths_df %>%
    arrange(chrom_index) %>%
    mutate(chr_offset = cumsum(as.numeric(lag(chr_len, default = 0))))
  
  # Merge chromosome lengths and offsets into asedf
  asedf_joined <- left_join(asedf, chrom_lengths_df[, c("contig", "chr_len", "chr_offset")], by = "contig")
  
  # Compute cumulative position
  asedf_joined <- asedf_joined %>%
    mutate(cumulative_pos = position + chr_offset)
  
  # Define significance threshold
  sig_threshold <- -log10(0.05)
  
  # Create the Manhattan plot
  p <- ggplot(data = asedf_joined, aes(x = cumulative_pos, y = -log10(pval))) +
    geom_point(aes(color = chrom_index %% 2 == 0), alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("skyblue", "navy")) +
    geom_hline(yintercept = sig_threshold, color = "red", linetype = "dashed") +
    labs(x = "Chromosome", y = "-log10(p-value)", title = "ASE Manhattan Plot") +
    theme_bw() +  # Use theme_bw() for white background
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    ylim(0, 10)  # Set the y-axis limit to a maximum of 10
  
  # Add chromosome labels
  axis_set <- chrom_lengths_df %>%
    mutate(center = chr_offset + chr_len / 2)
  
  p <- p + scale_x_continuous(
    breaks = axis_set$center,
    labels = axis_set$contig
  )
  
  # Annotate points above threshold with gene names
  asedf_joined_to_label <- asedf_joined %>%
    filter(-log10(pval) >= sig_threshold)
  
  p <- p + geom_text(
    data = asedf_joined_to_label,
    aes(x = cumulative_pos, y = -log10(pval), label = gene),
    size = 2, angle = 45, hjust = 0, nudge_x = 0, nudge_y = 0.1, check_overlap = TRUE
  )
  
  return(p)
}

run_ase_pipeline <- function(sample_id, reference_alleles_dir, ascat_counts_dir, results_dir, RNA_dir, WGS_dir, gatk_jar, java_cmd, java_home, gtf_file, min_depth = 3, filter_cutoff = 0.01, sig_threshold = -log10(0.05)) {
  
  # Step 1: Filter Allele Counts
  message("Filtering allele counts...")
  output_dir <- file.path(results_dir, sample_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  if(!dir.exists(file.path(output_dir, "allele_counts"))) {
    dir.create(file.path(output_dir, "allele_counts"), recursive = TRUE)
  }


  chromosomes <- c(1:22, "X")
  for (chr in chromosomes) {
    get_allele_counts_filtered(
      chr = chr,
      reference_alleles_dir = reference_alleles_dir,
      ascat_counts_dir = ascat_counts_dir,
      min_depth = min_depth,
      sample_id = sample_id,
      output_dir = output_dir
    )
  }
  
  # Step 2: Combine Loci Information
  message("Combining loci information...")
  combine_loci_nomatch(sample_id = sample_id, results_dir = results_dir)
  
  # Step 3: Run ASEReadCounter
  message("Running ASEReadCounter...")
  run_ASEReadCounter(
    sample_id = sample_id,
    results_dir = results_dir,
    RNA_dir = RNA_dir,
    WGS_dir = WGS_dir,
    gatk_jar = gatk_jar,
    java_cmd = java_cmd,
    java_home = java_home
  )
  
  # Step 4: Compute P-values
  message("Computing p-values...")
  asedf <- compute_pvals(
    sample_id = sample_id,
    results_dir = results_dir,
    filter_cutoff = filter_cutoff
  )
  
  # Step 5: Annotate Results
  message("Annotating results...")
  asedf_annotated <- annotate_ase_results(asedf, gtf_file = gtf_file)
  
  # Step 6: Plot Manhattan Plot
  message("Plotting Manhattan plot...")
  p <- plot_ase_manhattan(asedf_annotated, sig_threshold = sig_threshold)
  
  # Save plot
  plot_file <- file.path(results_dir, sample_id, paste0(sample_id, "_manhattan_plot.png"))
  ggsave(plot_file, plot = p, width = 12, height = 6, dpi = 300)
  
  message(paste("Manhattan plot saved to:", plot_file))
  
  return(asedf_annotated)
}

# Run the ASE analysis pipeline
# asedf_final <- run_ase_pipeline(
#   sample_id = sample_id,
#   reference_alleles_dir = reference_alleles_dir,
#   ascat_counts_dir = ascat_counts_dir,
#   results_dir = results_dir,
#   RNA_dir = RNA_dir,
#   WGS_dir = WGS_dir,
#   gatk_jar = gatk_jar,
#   java_cmd = java_cmd,
#   java_home = java_home,
#   gtf_file = gtf_file,
#   min_depth = 3,
#   filter_cutoff = 0.01,
#   sig_threshold = -log10(0.05)
# )

# Function to run ASE pipeline for all samples with DNA and RNA data
run_ase_pipeline_all_samples <- function(common_samples, reference_alleles_dir, ascat_counts_dir, results_dir, RNA_dir, WGS_dir, gatk_jar, java_cmd, java_home, gtf_file, min_depth = 3, filter_cutoff = 0.01, sig_threshold = -log10(0.05)) {
  for (sample_id in common_samples) {
    message(paste("Running ASE pipeline for sample:", sample_id))
    asedf <- run_ase_pipeline(
      sample_id = sample_id,
      reference_alleles_dir = reference_alleles_dir,
      ascat_counts_dir = ascat_counts_dir,
      results_dir = results_dir,
      RNA_dir = RNA_dir,
      WGS_dir = WGS_dir,
      gatk_jar = gatk_jar,
      java_cmd = java_cmd,
      java_home = java_home,
      gtf_file = gtf_file,
      min_depth = min_depth,
      filter_cutoff = filter_cutoff,
      sig_threshold = sig_threshold
    )
  }
}

# Checks for which manhattan plots are generated
finished_plots <- list.files(results_dir, full.names = TRUE, recursive = TRUE, pattern = "_manhattan_plot.png$")
finished_samples <- basename(dirname(finished_plots))

# Get list of samples that have not been processed
unfinished_samples <- common_samples[!common_samples %in% finished_samples]


# Run the ASE pipeline for all samples
run_ase_pipeline_all_samples(
  common_samples = unfinished_samples,
  reference_alleles_dir = reference_alleles_dir,
  ascat_counts_dir = ascat_counts_dir,
  results_dir = results_dir,
  RNA_dir = RNA_dir,
  WGS_dir = WGS_dir,
  gatk_jar = gatk_jar,
  java_cmd = java_cmd,
  java_home = java_home,
  gtf_file = gtf_file,
  min_depth = 3,
  filter_cutoff = 0.01,
  sig_threshold = -log10(0.05)
)
