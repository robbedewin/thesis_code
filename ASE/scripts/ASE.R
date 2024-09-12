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

get_allele_counts_filtered <- function(chr, reference_alleles_dir, sample_allele_counts_dir, mindepth = 3, sample_id, alias = "tumor") {
  # Conversion map
  conversion <- c("1" = "A", "2" = "C", "3" = "G", "4" = "T")
  
  # Construct file paths
  reference_alleles_file <- file.path(reference_alleles_dir, paste0("allele_T2T_chr", chr, ".txt"))
  sample_allele_counts_file <- file.path(sample_allele_counts_dir, paste0(sample_id, "/alleleFrequencies/", sample_id, "_", alias, "_alleleFrequencies_chr", chr, ".txt"))
  
  # Check if input files exist
  if (!file.exists(reference_alleles_file)) {
    warning(paste("Reference alleles file not found:", reference_alleles_file))
    return(NULL)
  }
  if (!file.exists(sample_allele_counts_file)) {
    warning(paste("Sample allele counts file not found:", sample_allele_counts_file))
    return(NULL)
  }
  
  # Output directory and file
  output_dir <- file.path("/staging/leuven/stg_00096/home/rdewin/ASE/results", sample_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  outfile_allelecounts <- file.path(output_dir, paste0(sample_id, "_", alias, "_filtered_allele_counts_chr", chr, ".txt"))
  
  # Read and process reference alleles
  reference_alleles <- tryCatch(
    {
      read_tsv(file = reference_alleles_file, col_names = TRUE, col_types = "iii")
    },
    error = function(e) {
      warning(paste("Error reading reference alleles file:", reference_alleles_file))
      return(NULL)
    }
  )
  
  if (is.null(reference_alleles)) return(NULL)
  
  colnames(reference_alleles) <- c("pos", "ref", "alt")
  
  # Convert numeric alleles to nucleotides
  reference_alleles$ref <- conversion[as.character(reference_alleles$ref)]
  reference_alleles$alt <- conversion[as.character(reference_alleles$alt)]
  
  # Handle potential NA values after conversion
  if (any(is.na(reference_alleles$ref)) | any(is.na(reference_alleles$alt))) {
    warning("NA values found after allele conversion. Check allele codes.")
    reference_alleles <- reference_alleles[!is.na(reference_alleles$ref) & !is.na(reference_alleles$alt), ]
  }
  
  # Read sample allele counts
  sample_allele_counts <- tryCatch(
    {
      read_tsv(file = sample_allele_counts_file, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
    },
    error = function(e) {
      warning(paste("Error reading sample allele counts file:", sample_allele_counts_file))
      return(NULL)
    }
  )
  
  if (is.null(sample_allele_counts)) return(NULL)
  
  # Combine chromosome information
  reference_alleles <- cbind(chr = sample_allele_counts$chr, reference_alleles)
  
  # Calculate counts for ref and alt alleles
  reference_alleles$count_ref <- with(reference_alleles, ifelse(ref == "A", count_A,
                                    ifelse(ref == "C", count_C, 
                                    ifelse(ref == "G", count_G, count_T))))
  
  reference_alleles$count_alt <- with(reference_alleles, ifelse(alt == "A", count_A,
                                    ifelse(alt == "C", count_C, 
                                    ifelse(alt == "G", count_G, count_T))))
  
  # Handle potential NA values after count calculation
  if (any(is.na(reference_alleles$count_ref)) | any(is.na(reference_alleles$count_alt))) {
    warning("NA values found in count_ref or count_alt. Check allele counts.")
    reference_alleles <- reference_alleles[!is.na(reference_alleles$count_ref) & !is.na(reference_alleles$count_alt), ]
  }
  
  # Filter based on minimum depth
  filtered_output <- subset(reference_alleles, count_ref >= mindepth & count_alt >= mindepth)
  
  if (nrow(filtered_output) == 0) {
    warning(paste("No loci passed the minimum depth filter for sample", sample_id, "chr", chr))
    return(NULL)
  }
  
  # Write to output file without column names
  tryCatch(
    {
      write_tsv(x = filtered_output, file = outfile_allelecounts, col_names = FALSE)
      message(paste("Filtered allele counts written to:", outfile_allelecounts))
    },
    error = function(e) {
      warning(paste("Error writing filtered allele counts to:", outfile_allelecounts))
    }
  )
  
  return(NULL)
}


# Define the parameters
chr <- "1"  # The chromosome number for which you want to get the allele counts
reference_alleles_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/ASCAT/ReferenceFiles"  # Directory where the reference allele files are stored
sample_allele_counts_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/ASCAT"  # Directory where sample-specific allele count files are stored
minimum_depth <- 5  # The minimum depth required for filtering the alleles
sample_id <- "P011"  # Unique identifier for the sample being analyzed
alias <- "tumor"  # Optional alias (e.g., "tumor" or "normal")

# Run the function
get_allele_counts_filtered(
  chr = chr,
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


combine_loci_nomatch_1 <- function(countsdir, sample_id, alias = "tumor") {
  # Define the sample-specific directory
  sample_dir <- file.path(countsdir, sample_id)
  
  # Check if sample directory exists
  if (!dir.exists(sample_dir)) {
    warning(paste("Sample directory not found:", sample_dir))
    return(NULL)
  }
  
  # List all matching files for the given sample_id and alias across all chromosomes
  allele_count_files <- list.files(path = sample_dir,
                                   pattern = paste0("^", sample_id, "_", alias, "_filtered_allele_counts_chr.*\\.txt$"),
                                   full.names = TRUE)
  
  if (length(allele_count_files) == 0) {
    warning(paste("No allele count files found for sample", sample_id, "with alias", alias))
    return(NULL)
  }
  
  # Combine the allele count files into a single data frame
  combined_allele_counts <- do.call(rbind, lapply(X = allele_count_files, FUN = function(file) {
    tryCatch(
      {
        read_tsv(file, col_names = c("chromosome", "position", "reference_allele", "alternate_allele", "reference_allele_count", "alternate_allele_count"), col_types = "ciccii")
      },
      error = function(e) {
        warning(paste("Error reading file:", file))
        return(NULL)
      }
    )
  }))
  
  # Remove NULL entries resulting from failed reads
  combined_allele_counts <- combined_allele_counts[!sapply(combined_allele_counts, is.null)]
  
  if (nrow(combined_allele_counts) == 0) {
    warning("No valid allele count data to combine.")
    return(NULL)
  }
  
  # Write the combined allele counts to a text file
  output_file_txt <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  tryCatch(
    {
      write_tsv(x = combined_allele_counts, file = output_file_txt, col_names = TRUE)
      message(paste("Combined allele counts written to:", output_file_txt))
    },
    error = function(e) {
      warning(paste("Error writing combined allele counts to:", output_file_txt))
    }
  )
  
  # Write the combined allele counts to a VCF file
  output_file_vcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))
  
  # Determine if chromosome names in the file start with "chr"
  sample_size <- min(100, nrow(combined_allele_counts))
  chromosome_has_prefix <- any(grepl(pattern = "^chr", x = combined_allele_counts$chromosome[sample(x = 1:nrow(combined_allele_counts), size = sample_size, replace = TRUE)]))
  
  # Create VRanges based on chromosome naming
  allelecounts_vr <- tryCatch(
    {
      if (chromosome_has_prefix) {
        VRanges(seqnames = combined_allele_counts$chromosome,
                ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                ref = combined_allele_counts$reference_allele,
                alt = combined_allele_counts$alternate_allele,
                seqinfo = seqinfo(bsgenome_hs1))
      } else {
        VRanges(seqnames = paste0("chr", combined_allele_counts$chromosome),
                ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                ref = combined_allele_counts$reference_allele,
                alt = combined_allele_counts$alternate_allele,
                seqinfo = seqinfo(bsgenome_T2T))
      }
    },
    error = function(e) {
      warning("Error creating VRanges object.")
      return(NULL)
    }
  )
  
  if (is.null(allelecounts_vr)) return(NULL)
  
  # Check for NA seqnames
  if (any(is.na(seqnames(allelecounts_vr)))) {
    warning("NA values found in seqnames. Check chromosome naming conventions.")
    return(NULL)
  }
  
  # Create matrices for genotype fields (DP = depth, AD = allele depth)
  depth_matrix <- as.integer(combined_allele_counts$reference_allele_count + combined_allele_counts$alternate_allele_count)
  allele_depth_matrix <- cbind(as.integer(combined_allele_counts$reference_allele_count), as.integer(combined_allele_counts$alternate_allele_count))
  
  # Ensure matrices have correct dimensions
  if (length(depth_matrix) != length(allelecounts_vr) | nrow(allele_depth_matrix) != length(allelecounts_vr)) {
    warning("Mismatch in dimensions between allele counts and VRanges object.")
    return(NULL)
  }
  
  # Add genotype fields to VRanges object
  geno(allelecounts_vr) <- SimpleList(
    DP = matrix(depth_matrix, ncol = 1, dimnames = list(NULL, "DP")),
    AD = allele_depth_matrix
  )
  
  # Assign column names for AD
  colnames(geno(allelecounts_vr)$AD) <- c(paste0(sample_id, "_ref"), paste0(sample_id, "_alt"))
  
  # Sort the VRanges object and assign the sample name
  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  
  # Write the VRanges object to a VCF file
  tryCatch(
    {
      writeVcf(obj = allelecounts_vr, filename = output_file_vcf, index = TRUE)
      message(paste("VCF file written to:", output_file_vcf))
    },
    error = function(e) {
      warning(paste("Error writing VCF file:", output_file_vcf))
    }
  )
  
  return(NULL)
}




process_sample_ase <- function(sample_id, alias = "tumor", chromosomes = c(1:22, "X"),
                               reference_alleles_dir, sample_allele_counts_dir, countsdir) {
  # Step 1: Filter Allele Counts per Chromosome
  message(paste("Starting ASE processing for sample:", sample_id, "Alias:", alias))
  for (chr in chromosomes) {
    message(paste("Processing chromosome:", chr))
    get_allele_counts_filtered(
      chr = chr,
      reference_alleles_dir = reference_alleles_dir,
      sample_allele_counts_dir = sample_allele_counts_dir,
      mindepth = 5,  # Adjust based on your requirements
      sample_id = sample_id,
      alias = alias
    )
  }
  
  # Step 2: Combine Filtered Allele Counts
  message("Combining filtered allele counts across chromosomes.")
  combine_loci_nomatch_1(
    countsdir = countsdir,
    sample_id = sample_id,
    alias = alias
  )
  
  # Proceed to Step 3: Run ASEReadCounter (Function not yet provided)
  # Example:
  # run_ase_read_counter(sample_id)
  
  # Further Steps: Compute p-values, Annotate, Plot (Functions to be implemented)
  
  message(paste("Completed ASE processing for sample:", sample_id))
}
