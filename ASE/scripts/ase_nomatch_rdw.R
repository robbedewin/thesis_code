# Set library path
.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

# Install necessary packages if not already installed
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hs1")
# install.packages("readr")
# install.packages("rslurm")
# install.packages("readr")
# BiocManager::install("VariantAnnotation")
# BiocManager::install("biomaRt")

# Load necessary libraries
library(readr)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(BSgenome.Hsapiens.UCSC.hs1)
library(biomaRt)
library(ggplot2)
library(parallel)
# library(VGAM)

# Load BSgenome data
bsgenome_T2T <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
bsgenome_hs1 <- BSgenome.Hsapiens.UCSC.hs1

# Source utility scripts
source(file = "/staging/leuven/stg_00096/home/rdewin/WGS/rules/scripts/1000Genomes_getAllelecounts.R")
source(file = "/staging/leuven/stg_00096/home/rdewin/WGS/rules/scripts/utils.R")

# Define paths to required software and reference files
ALLELECOUNTER <- "/staging/leuven/stg_00096/software/alleleCount/4.0.0-GCCcore-6.4.0/bin/alleleCounter"
ALLELESDIR <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/ReferenceFiles/1000G/"
JAVA <- "/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-2.el8.x86_64/jre/bin/java"
ALLELESDIR <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/ReferenceFiles/"


# Function to get allelecount for ref and alt and filter based on minimum depth


# Function to get allele count for reference and alternate alleles and filter based on minimum depth

get_allele_counts_filtered <- function(chr, reference_alleles_dir, sample_allele_counts_dir, mindepth = 3, sample_id, alias = "tumor") {
  # Define the conversion map for numeric to nucleotide conversion
  conversion <- c("1" = "A", "2" = "C", "3" = "G", "4" = "T")
  
  # Construct file paths for reference alleles and sample allele counts
  reference_alleles_file <- file.path(reference_alleles_dir, paste0("allele_T2T_chr", chr, ".txt"))
  sample_allele_counts_file <- file.path(sample_allele_counts_dir, paste0(sample_id, "/alleleFrequencies/", sample_id, "_", alias, "_alleleFrequencies_chr", chr, ".txt"))
    
  # Define the output directory for filtered allele counts
  output_dir <- file.path("/staging/leuven/stg_00096/home/rdewin/ASE/results", sample_id)
  
  # Ensure the output directory exists, create it if necessary
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define the output file path for filtered allele counts
  outfile_allelecounts <- file.path(output_dir, paste0(sample_id, "_", alias, "_filtered_allele_counts_chr", chr, ".txt"))
  
  
  # Read the reference alleles data with original column names
  reference_alleles <- read_tsv(file = reference_alleles_file, col_names = TRUE, col_types = "iii")
  
  # Rename the columns to c("pos", "ref", "alt")
  colnames(reference_alleles) <- c("pos", "ref", "alt")
  
  # Convert numeric alleles to nucleotide bases
  reference_alleles$ref <- conversion[as.character(reference_alleles$ref)]
  reference_alleles$alt <- conversion[as.character(reference_alleles$alt)]
  
  # Read the sample allele counts data with specified column types
  sample_allele_counts <- read_tsv(file = sample_allele_counts_file, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
  
  # Combine the chromosome name with the reference alleles data
  reference_alleles <- cbind(chr = sample_allele_counts$chr, reference_alleles)
  
  # Calculate reference and alternate allele counts
  reference_alleles$count_ref <- ifelse(reference_alleles$ref == "A", sample_allele_counts$count_A,
                                        ifelse(reference_alleles$ref == "C", sample_allele_counts$count_C, 
                                               ifelse(reference_alleles$ref == "G", sample_allele_counts$count_G, sample_allele_counts$count_T)))
  reference_alleles$count_alt <- ifelse(reference_alleles$alt == "A", sample_allele_counts$count_A,
                                        ifelse(reference_alleles$alt == "C", sample_allele_counts$count_C, 
                                               ifelse(reference_alleles$alt == "G", sample_allele_counts$count_G, sample_allele_counts$count_T)))
  
  # Filter the data based on the minimum depth for both reference and alternate allele counts
  filtered_output <- reference_alleles[reference_alleles$count_ref >= mindepth & reference_alleles$count_alt >= mindepth, ]
  
  # Write the filtered data to an output file
  write_tsv(x = filtered_output, file = outfile_allelecounts, col_names = F)

  
  return(NULL)
}


# Define the parameters
chromosome <- "1"  # The chromosome number for which you want to get the allele counts
reference_alleles_directory <- "/staging/leuven/stg_00096/home/rdewin/ASE/ASCAT/ReferenceFiles"  # Directory where the reference allele files are stored
sample_allele_counts_directory <- "/staging/leuven/stg_00096/home/rdewin/ASE/ASCAT"  # Directory where sample-specific allele count files are stored
minimum_depth <- 5  # The minimum depth required for filtering the alleles
sample_id <- "P011"  # Unique identifier for the sample being analyzed
alias <- "tumor"  # Optional alias (e.g., "tumor" or "normal")

# Run the function
get_allele_counts_filtered(
  chr = chromosome,
  reference_alleles_dir = reference_alleles_directory,
  sample_allele_counts_dir = sample_allele_counts_directory,
  mindepth = minimum_depth,
  sample_id = sample_id,
  alias = alias
)

chromosomes <- c(1:22, "X")
for (chr in chromosomes) {
  # Call the get_allele_counts_filtered function for each chromosome
  get_allele_counts_filtered(chr = chr, 
                             reference_alleles_dir = reference_alleles_directory,
                             sample_allele_counts_dir = sample_allele_counts_directory,
                             mindepth = minimum_depth,
                             sample_id = sample_id,
                             alias = alias)
  
  # Optionally print a message to track progress
  print(paste("Processed chromosome", chr))
}


# This function combines allele counts and creates a VCF file using VRanges, handling chromosome name formats.
combine_loci_nomatch_1 <- function(countsdir, sample_id, alias = "tumor") {
  # List all matching files for the given sample_id and alias across all chromosomes
  allele_count_files <- list.files(path = file.path(countsdir, sample_id),
                                   pattern = paste0(sample_id, "_", alias, "_filtered_allele_counts_chr", ".*\\.txt$"),
                                   full.names = TRUE)

  # Combine the allele count files into a single data frame
  combined_allele_counts <- do.call(rbind, lapply(X = allele_count_files, FUN = read_tsv, col_names = c("chromosome", "position", "reference_allele", "alternate_allele", "reference_allele_count", "alternate_allele_count"), col_types = "ciccii"))
  
  # Write the combined allele counts to a text file
  output_file_txt <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  write_tsv(x = combined_allele_counts, file = output_file_txt, col_names = TRUE)
  
  # Write the combined allele counts to a VCF file
  output_file_vcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))

  # Determine if chromosome names in the file start with "chr" and adjust if necessary
  chromosome_has_prefix <- any(grepl(pattern = "^chr", x = combined_allele_counts$chromosome[sample(x = 1:nrow(combined_allele_counts), size = 100, replace = TRUE)]))
  
  # Create VRanges based on whether chromosome names have "chr" prefix
  if (chromosome_has_prefix) {
    allelecounts_vr <- VRanges(seqnames = combined_allele_counts$chromosome,
                               ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                               ref = combined_allele_counts$reference_allele, alt = combined_allele_counts$alternate_allele,
                               seqinfo = seqinfo(bsgenome_hs1))
  } else {
    allelecounts_vr <- VRanges(seqnames = paste0("chr", combined_allele_counts$chromosome),
                               ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                               ref = combined_allele_counts$reference_allele, alt = combined_allele_counts$alternate_allele,
                               seqinfo = seqinfo(bsgenome_T2T))
  }

  # Create matrices for genotype fields (DP = depth, AD = allele depth)
  depth_matrix <- matrix(as.integer(combined_allele_counts$reference_allele_count + combined_allele_counts$alternate_allele_count), nrow = length(allelecounts_vr), ncol = 1)
  allele_depth_matrix <- matrix(c(combined_allele_counts$reference_allele_count, combined_allele_counts$alternate_allele_count), nrow = length(allelecounts_vr), ncol = 2)

  # Set column names for the matrices
  colnames(depth_matrix) <- sample_id
  colnames(allele_depth_matrix) <- c(paste0(sample_id, "_ref"), paste0(sample_id, "_alt"))

  # Add genotype fields to VRanges object
  geno(allelecounts_vr) <- SimpleList(
    DP = depth_matrix,
    AD = allele_depth_matrix
  )

  # Sort the VRanges object and assign the sample name
  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  
  # Write the VRanges object to a VCF file
  writeVcf(obj = allelecounts_vr, filename = output_file_vcf, index = TRUE)
  
  return(NULL)
}

countsdir = "/staging/leuven/stg_00096/home/rdewin/ASE/results"
sample_id = "P011"
alias = "tumor"

# Example usage for a tumor sample
combine_loci_nomatch_1(
  countsdir = "/staging/leuven/stg_00096/home/rdewin/ASE/results",
  sample_id = "P011",
  alias = "tumor"
)



# This function combines allele counts and creates a VCF file using GRanges and direct VCF creation, ensuring chromosome name compatibility.
combine_loci_nomatch_2 <- function(countsdir, sample_id, alias="tumor") {
  # List files in the counts directory matching the pattern
  allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_nomatch_chr", full.names = TRUE)

  # Combine the allele counts files
  allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))
  
  # Write the combined allele counts to a file
  outfile <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  write_tsv(x = allelecounts, path = outfile, col_names = TRUE)
  
  # Write the combined allele counts to a VCF file
  locivcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))

  # Determine if the chromosome names start with "chr"
  if (any(grepl(pattern = "^chr", x = allelecounts$chr[sample(x = 1:nrow(allelecounts), size = 100, replace = TRUE)]))) {
    seqinfo <- seqinfo(bsgenome_hs1)
  } else {
    seqinfo <- seqinfo(bsgenome_T2T)
    allelecounts$chr <- paste0("chr", allelecounts$chr)
  }

  # Create GRanges object
  gr <- GRanges(seqnames = allelecounts$chr,
                ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos),
                seqinfo = seqinfo)
  
  # Create colData DataFrame
  col_data <- DataFrame(Samples = sample_id, row.names = sample_id)
  
  # Create the DP and AD matrices with the correct dimensions
  dp_matrix <- matrix(as.integer(allelecounts$count_ref + allelecounts$count_alt), nrow = length(gr), ncol = 1)
  ad_matrix <- matrix(c(allelecounts$count_ref, allelecounts$count_alt), nrow = length(gr), ncol = 2)

  # Ensure the column names of the matrices are set
  colnames(dp_matrix) <- sample_id
  colnames(ad_matrix) <- c(paste0(sample_id, "_R"), paste0(sample_id, "_A"))

  # Create VCF object
  vcf <- VCF(rowRanges = gr,
             colData = col_data,
             geno = SimpleList(
               DP = dp_matrix,
               AD = ad_matrix
             ))

  # Write VCF to file
  writeVcf(vcf, locivcf)
  
  return(NULL)
}

# Execute the functions
sample_id <- "P011"
countsdir <- "/staging/leuven/stg_00096/home/rdewin/ASE/"
combine_loci_nomatch_2(countsdir = countsdir, sample_id = sample_id)



combine_loci_nomatch_3 <- function(countsdir, sample_id) {
  
  # List files in the counts directory matching the pattern
  allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_nomatch_chr", full.names = TRUE)

  # Combine the allele counts files
  allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))

  # Write the combined allele counts to a file
  outfile <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  write_tsv(x = allelecounts, path = outfile, col_names = TRUE)
  
  # Write the combined allele counts to a VCF file
  locivcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))

  # Determine if the chromosome names start with "chr"
  if (any(grepl(pattern = "^chr", x = allelecounts$chr[sample(x = 1:nrow(allelecounts), size = 100, replace = TRUE)]))) {
    allelecounts_vr <- VRanges(seqnames = allelecounts$chr,
                               ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos),
                               ref = allelecounts$ref, alt = allelecounts$alt,
                               seqinfo = seqinfo(bsgenome_hs1))
  } else {
    allelecounts_vr <- VRanges(seqnames = paste0("chr", allelecounts$chr),
                               ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos),
                               ref = allelecounts$ref, alt = allelecounts$alt,
                               seqinfo = seqinfo(bsgenome_T2T))
  }
  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  writeVcf(obj = allelecounts_vr, filename = locivcf, index = TRUE)
  
  return(NULL)
}


# Execute the functions
sample_id <- "P011"
countsdir <- "/staging/leuven/stg_00096/home/rdewin/ASE/"
#get_alleles_chr_nomatch(chr = "1", allelesdir = ALLELESDIR, countsdir = countsdir, sample_id = sample_id)
combine_loci_nomatch_1(countsdir = countsdir, sample_id = sample_id)
combine_loci_nomatch_2(countsdir = countsdir, sample_id = sample_id)
combine_loci_nomatch_3(countsdir = countsdir, sample_id = sample_id)

# Use ASEReadCounter to get allele counts on the RNA
ASEReadCount <- function(hetSNPvcf, bamfile, refgenome, outfile, minBaseQ = 20, minMapQ = 35) {
  
  # Path to GATK executable
  gatk.exe = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/gatk"
  
  # Path to Java executable in Conda environment
  java_cmd <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/java"
  
  # Path to GATK jar file
  gatk_jar <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar"
  
  # Set JAVA_HOME to use the correct Java version
  Sys.setenv(JAVA_HOME = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS")
  
  # Construct the command
  cmd <- paste0(java_cmd, " -Xmx12G -jar ", gatk_jar, " ASEReadCounter",
                " -I ", bamfile,
                " -V ", hetSNPvcf,
                " -R ", refgenome,
                " -O ", outfile,
                " --min-depth 1",
                " --min-mapping-quality ", minMapQ,
                " --min-base-quality ", minBaseQ,
                " --lenient")
  
  cat("Executing command:\n", cmd, "\n")
  system(cmd, wait = TRUE)
}



# Define file paths and GATK executable
hetSNPvcf <- "/staging/leuven/stg_00096/home/rdewin/ASE/P011_hetSNPs_nomatch.vcf.bgz"
bamfile <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/star/P011/P011_Aligned.sortedByCoord.withRG.bam"
refgenome <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa"
outfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/P011_asereadcounts_nomatch.rtable"
gatk <- "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/gatk"

# Execute the function
ASEReadCount(hetSNPvcf = hetSNPvcf, bamfile = bamfile, refgenome = refgenome, outfile = outfile)
