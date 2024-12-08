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




# Source the beta-binomial test function
source(file = "//staging/leuven/stg_00096/home/rdewin/ASE/scripts/utils.R") 

# Source the ASE pipeline functions
source(file = "/staging/leuven/stg_00096/home/rdewin/ASE/scripts/ase_analysis_functions.R")

# Load the libraries
load_libraries()

# Set variables defined in the analysis_function script
vars <- set_ase_variables()

# Assign variables to the global environment
list2env(vars, .GlobalEnv)

# Get common samples
common_samples <- get_common_samples(RNA_dir, ascat_counts_dir)


# Ensure output directories exist
for (sample_id in common_samples){
  create_output_dirs(results_dir, sample_id)
}

# Function to filter allele counts
for (sample_id in common_samples) {
  for (chr in chromosomes) {
    get_allele_counts_filtered(
      chr = chr,
      reference_alleles_dir = reference_alleles_dir,
      ascat_counts_dir = ascat_counts_dir,
      min_depth = 3,
      sample_id = sample_id,
      output_dir = file.path(results_dir, sample_id)
    )
  }
}

# Function to combine loci information
for (sample_id in common_samples) {
  combine_loci_nomatch(sample_id = sample_id, results_dir = results_dir)
}

# Check if ASEReadCounter already ran
aserun <- list.files(results_dir, full.names = TRUE, recursive = TRUE, pattern = "asereadcounts_nomatch.tsv$")
aserun_samples <- basename(dirname(aserun))
samples_not_ASEReadCounter <- common_samples[!common_samples %in% aserun_samples]

# Function to run ASEReadCounter
for (sample_id in samples_not_ASEReadCounter) {
  run_ASEReadCounter(
    sample_id = sample_id,
    results_dir = results_dir,
    RNA_dir = RNA_dir,
    ref_genome = ref_genome,
    gatk_jar = gatk_jar,
    java_cmd = java_cmd,
    java_home = java_home
  )
}

# Function to compute p-values
for (sample_id in common_samples) {
  message(paste("Computing p-values for sample:", sample_id))
  asedf <- compute_pvals(
    sample_id = sample_id,
    results_dir = results_dir,
    filter_cutoff = 0.01
  )
}

# Alternative function to compute p-values
for (sample_id in common_samples) {
  asedf <- compute_pvals_alternative(
    sample_id = sample_id,
    results_dir = results_dir,
    filter_cutoff = 0.01
  )
}

# Function to plot QC metrics
for (sample_id in common_samples) {
  plot_qc_plots(results_dir, sample_id)
}     

# Function to annotate ASE results
for (sample_id in common_samples) {
  annotate_ase_results(sample_id, results_dir, gtf_file)
}

# Function to plot Manhattan plot
for (sample_id in common_samples) {
  message(paste("Plotting Manhattan plot for sample:", sample_id))
  p <- plot_ase_manhattan(sample_id, results_dir, sig_threshold = -log10(0.05))
  
  # Save plot
  plot_file <- file.path(results_dir, sample_id, paste0(sample_id, "_manhattan_plot_new.png"))
  ggsave(plot_file, plot = p, width = 12, height = 6, dpi = 300)
  
  message(paste("Manhattan plot saved to:", plot_file))
}




run_ase_pipeline <- function(sample_id, chromosomes, reference_alleles_dir, ascat_counts_dir, results_dir, RNA_dir, ref_genome, gatk_jar, java_cmd, java_home, gtf_file, min_depth = 3, filter_cutoff = 0.01, sig_threshold = -log10(0.05)) {
  
  # Step 1: Filter Allele Counts
  message("Filtering allele counts...")
  
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
    ref_genome = ref_genome,
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



# Function to run ASE pipeline for all samples with DNA and RNA data
run_ase_pipeline_all_samples <- function(common_samples, chromosomes, reference_alleles_dir, ascat_counts_dir, results_dir, RNA_dir, ref_genome, gatk_jar, java_cmd, java_home, gtf_file, min_depth = 3, filter_cutoff = 0.01, sig_threshold = -log10(0.05)) {
  for (sample_id in common_samples) {
    message(paste("Running ASE pipeline for sample:", sample_id))
    asedf <- run_ase_pipeline(
      sample_id = sample_id,
      chromosomes = chromosomes,
      reference_alleles_dir = reference_alleles_dir,
      ascat_counts_dir = ascat_counts_dir,
      results_dir = results_dir,
      RNA_dir = RNA_dir,
      ref_genome = ref_genome,
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
  ref_genome = ref_genome,
  gatk_jar = gatk_jar,
  java_cmd = java_cmd,
  java_home = java_home,
  gtf_file = gtf_file,
  min_depth = 3,
  filter_cutoff = 0.01,
  sig_threshold = -log10(0.05)
)
