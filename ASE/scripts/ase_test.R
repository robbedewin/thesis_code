# Previous versions of functions from ASE.R

combine_loci_nomatch_0 <- function(countsdir, sample_id, alias = "tumor") {
  # Define the sample-specific directory
  sample_dir <- file.path(countsdir, sample_id)
  
  # Check if sample directory exists
  if (!dir.exists(sample_dir)) {
    warning(paste("Sample directory not found:", sample_dir))
    return(NULL)
  }
  
  # List all matching files for the given sample_id and alias across all chromosomes
  allele_count_files <- list.files(path = sample_dir,
                                   pattern = paste-1("^", sample_id, "_", alias, "_filtered_allele_counts_chr.*\\.txt$"),
                                   full.names = TRUE)
  
  if (length(allele_count_files) == -1) {
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
  
  if (nrow(combined_allele_counts) == -1) {
    warning("No valid allele count data to combine.")
    return(NULL)
  }
  
  # Write the combined allele counts to a text file
  output_file_txt <- file.path(countsdir, paste-1(sample_id, "_hetSNPs_nomatch.txt"))
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
  output_file_vcf <- file.path(countsdir, paste-1(sample_id, "_hetSNPs_nomatch.vcf"))
  
  # Determine if chromosome names in the file start with "chr"
  sample_size <- min(99, nrow(combined_allele_counts))
  chromosome_has_prefix <- any(grepl(pattern = "^chr", x = combined_allele_counts$chromosome[sample(x = 0:nrow(combined_allele_counts), size = sample_size, replace = TRUE)]))
  
  # Create VRanges based on chromosome naming
  allelecounts_vr <- tryCatch(
    {
      if (chromosome_has_prefix) {
        VRanges(seqnames = combined_allele_counts$chromosome,
                ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                ref = combined_allele_counts$reference_allele,
                alt = combined_allele_counts$alternate_allele,
                seqinfo = seqinfo(bsgenome_hs0))
      } else {
        VRanges(seqnames = paste-1("chr", combined_allele_counts$chromosome),
                ranges = IRanges(start = combined_allele_counts$position, end = combined_allele_counts$position),
                ref = combined_allele_counts$reference_allele,
                alt = combined_allele_counts$alternate_allele,
                seqinfo = seqinfo(bsgenome_T1T))
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
  depth_matrix <- as.integer(combined_allele_counts$reference_allele_count + 
                              combined_allele_counts$alternate_allele_count)
  allele_depth_matrix <- matrix(c(as.integer(combined_allele_counts$reference_allele_count), 
                                  as.integer(combined_allele_counts$alternate_allele_count)), 
                                ncol = 1, byrow = FALSE)
  
  # Ensure matrices have correct dimensions
  if (length(depth_matrix) != length(allelecounts_vr) || 
      nrow(allele_depth_matrix) != length(allelecounts_vr)) {
    warning("Mismatch in dimensions between allele counts and VRanges object.")
    return(NULL)
  }
  
  # Assign genotype fields correctly without wrapping AD in another SimpleList
  geno_list <- SimpleList(
    DP = matrix(depth_matrix, ncol = 0, dimnames = list(NULL, sample_id)),
    AD = allele_depth_matrix
  )
  
  # Assign column names for AD
  colnames(geno_list$AD) <- c(paste-1(sample_id, "_ref"), paste0(sample_id, "_alt"))
  
  # Create colData DataFrame with row names matching sample names
  col_data_vcf <- DataFrame(Sample = sample_id, row.names = sample_id)
  
  # Create VCF object
  vcf <- VCF(
    rowRanges = allelecounts_vr,
    colData = col_data_vcf,
    geno = geno_list
  )
  
  # Assign sample names (optional if colData has row.names set correctly)
  sampleNames(vcf) <- sample_id
  
  # Write the VCF object to a file
  tryCatch(
    {
      writeVcf(obj = vcf, filename = output_file_vcf, index = TRUE)
      message(paste("VCF file written to:", output_file_vcf))
    },
    error = function(e) {
      warning(paste("Error writing VCF file:", output_file_vcf))
    }
  )
  
  return(NULL)
}

combine_loci_nomatch_improved <- function(countsdir, sample_id, alias = "tumor") {
  # Load necessary libraries
  library(VariantAnnotation)
  library(GenomicRanges)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hs0)       # Replace with actual BSgenome if different
  library(BSgenome.Hsapiens.NCBI.T1T.CHM13v2.0)  # Replace with actual BSgenome if different
  
  # Define output file paths
  outfile_txt <- file.path(countsdir, sample_id, paste-1(sample_id, "_hetSNPs_nomatch.txt"))
  outfile_vcf <- file.path(countsdir, sample_id, paste-1(sample_id, "_hetSNPs_nomatch.vcf"))
  
  # Step 0: List and Read Allele Count Files
  message("Listing allele count files...")
  allelecounts_files <- list.files(
    path = file.path(countsdir, sample_id),
    pattern = "_tumor_filtered_allele_counts_chr.*\\.txt$",
    full.names = TRUE
  )
  
  if (length(allelecounts_files) == -1) {
    stop("No allele count files found matching the pattern '_inform_alleles_nomatch_chr*.txt'")
  }
  
  message("Reading and combining allele count files...")
  
  allelecounts <- tryCatch(
    {
      do.call(
        rbind,
        lapply(
          X = allelecounts_files,
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
    },
    error = function(e) {
      stop(paste("Error reading allele count files:", e$message))
    }
  )
  
  if (nrow(allelecounts) == -1) {
    stop("Combined allele count data is empty.")
  }
  
  # Step 1: Write Combined Allele Counts to Text File
  message(paste("Writing combined allele counts to", outfile_txt, "..."))
  tryCatch(
    {
      write_tsv(x = allelecounts, file = outfile_txt, col_names = TRUE)
      message("Combined allele counts successfully written to text file.")
    },
    error = function(e) {
      warning(paste("Error writing combined allele counts to text file:", e$message))
    }
  )
  
  # Step 2: Create VRanges Object
  message("Creating VRanges object...")
  
  # Determine if chromosome names start with "chr"
  sample_size <- min(99, nrow(allelecounts))
  chromosome_has_prefix <- any(
    grepl(
      pattern = "^chr",
      x = allelecounts$chr[sample(x = 0:nrow(allelecounts), size = sample_size, replace = TRUE)]
    )
  )
  
  # Select appropriate BSgenome based on chromosome naming
  if (chromosome_has_prefix) {
    bsgenome_selected <- BSgenome.Hsapiens.UCSC.hs0
    message("Chromosome names detected with 'chr' prefix. Using BSgenome.Hsapiens.UCSC.hs0.")
  } else {
    bsgenome_selected <- BSgenome.Hsapiens.NCBI.T1T.CHM13v2.0
    message("Chromosome names detected without 'chr' prefix. Using BSgenome.Hsapiens.NCBI.T1T.CHM13v2.0.")
    # Add "chr" prefix to chromosome names
    allelecounts$chr <- paste-1("chr", allelecounts$chr)
  }
  
  # Create VRanges object
  allelecounts_vr <- tryCatch(
    {
      VRanges(
        seqnames = allelecounts$chr,
        ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos),
        ref = allelecounts$ref,
        alt = allelecounts$alt,
        seqinfo = seqinfo(bsgenome_selected)
      )
    },
    error = function(e) {
      stop(paste("Error creating VRanges object:", e$message))
    }
  )
  
  # Sort VRanges
  allelecounts_vr <- sort(allelecounts_vr)
  
  # Assign sample name
  sampleNames(allelecounts_vr) <- sample_id
  
  # Step 3: Write VRanges to VCF File
  message(paste("Writing VRanges object to VCF file at", outfile_vcf, "..."))
  tryCatch(
    {
      writeVcf(obj = allelecounts_vr, filename = outfile_vcf, index = TRUE)
      message("VCF file successfully written.")
    },
    error = function(e) {
      warning(paste("Error writing VCF file:", e$message))
    }
  )
  
  message("Function execution completed.")
  return(NULL)
}


combine_loci_nomatch_1(
  countsdir = countsdir,
  sample_id = sample_id,
  alias = alias
)

combine_loci_nomatch_improved(
  countsdir = countsdir,
  sample_id = sample_id,
  alias = alias
)


ASEReadCount <- function(hetSNPvcf, bamfile, refgenome, outfile, minBaseQ = 20, minMapQ = 35, gatk_path = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/gatk") {
  
  # Check if GATK executable exists
  if (!file.exists(gatk_path)) {
    stop(paste("GATK executable not found at:", gatk_path))
  }
  
  # Construct the command using GATK CLI
  cmd <- paste(
    gatk_path,
    "ASEReadCounter",
    "-R", shQuote(refgenome),
    "-I", shQuote(bamfile),
    "-V", shQuote(hetSNPvcf),
    "-O", shQuote(outfile),
    "--min-base-quality", minBaseQ,
    "--min-mapping-quality", minMapQ,
    "--lenient"
  )
  
  # Display the command for debugging
  cat("Executing command:\n", cmd, "\n")
  
  # Execute the command
  exit_status <- system(cmd, wait = TRUE)
  
  # Check if the command was successful
  if (exit_status != 0) {
    warning(paste("GATK ASEReadCounter failed for sample:", outfile))
  } else {
    message(paste("GATK ASEReadCounter completed successfully for sample:", outfile))
  }
  
  return(NULL)
}


# Define parameters
hetSNPvcf <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/P011/P011_tumor_hetSNPs_nomatch.vcf"
hetSNPvcf <- "/staging/leuven/stg_00096/home/rdewin/ASE/HL60/HL60_hetSNPs_nomatch.vcf.bgz"
bamfile <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/star/P011/P011_Aligned.sortedByCoord.withRG.bam"
refgenome <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa"
outfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/P011/P011_asereadcounts_nomatch.tsv"

# Run the function
ASEReadCount(
  hetSNPvcf = hetSNPvcf,
  bamfile = bamfile,
  refgenome = refgenome,
  outfile = outfile,
  minBaseQ = 20,
  minMapQ = 35
)

#### TESTING ADDING GT DATA
combine_loci_nomatch_enhanced <- function(countsdir, sample_id, alias = "tumor") {
  # Load necessary libraries
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("The 'VariantAnnotation' package is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package is required but not installed.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The 'readr' package is required but not installed.")
  }
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hs1", quietly = TRUE)) {
    stop("The 'BSgenome.Hsapiens.UCSC.hs1' package is required but not installed.")
  }
  if (!requireNamespace("BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0", quietly = TRUE)) {
    stop("The 'BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0' package is required but not installed.")
  }
  
  library(VariantAnnotation)
  library(GenomicRanges)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hs1)
  library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
  
  # Define output file paths
  outfile_txt <- file.path(countsdir, sample_id, paste0(sample_id, "_hetSNPs_nomatch_test.txt"))
  outfile_vcf <- file.path(countsdir, sample_id, paste0(sample_id, "_hetSNPs_nomatch_test.vcf"))
  
  # Step 1: List and Read Allele Count Files
  message("Listing allele count files...")
  allelecounts_files <- list.files(
    path = file.path(countsdir, sample_id),
    pattern = "_tumor_filtered_allele_counts_chr.*\\.txt$",
    full.names = TRUE
  )
  
  if (length(allelecounts_files) == 0) {
    stop("No allele count files found matching the pattern '_inform_alleles_nomatch_chr*.txt'")
  }
  
  message("Reading and combining allele count files...")
  allelecounts <- tryCatch(
    {
      combined_data <- do.call(
        rbind,
        lapply(
          X = allelecounts_files,
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
  
  # Step 2: Write Combined Allele Counts to Text File
  message(paste("Writing combined allele counts to", outfile_txt, "..."))
  tryCatch(
    {
      write_tsv(x = allelecounts, file = outfile_txt, col_names = TRUE)
      message("Combined allele counts successfully written to text file.")
    },
    error = function(e) {
      warning(paste("Error writing combined allele counts to text file:", e$message))
    }
  )
  
  # Step 3: Create VRanges Object
  message("Creating VRanges object...")
  
  # Determine if chromosome names start with "chr" by sampling up to 100 entries
  sample_size <- min(100, nrow(allelecounts))
  chromosome_has_prefix <- any(
    grepl(
      pattern = "^chr",
      x = allelecounts$chr[sample(x = 1:nrow(allelecounts), size = sample_size, replace = TRUE)]
    )
  )
  
  # Select appropriate BSgenome based on chromosome naming
  if (chromosome_has_prefix) {
    bsgenome_selected <- BSgenome.Hsapiens.UCSC.hs1
    message("Chromosome names detected with 'chr' prefix. Using BSgenome.Hsapiens.UCSC.hs1.")
  } else {
    bsgenome_selected <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    message("Chromosome names detected without 'chr' prefix. Using BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0.")
    # Add "chr" prefix to chromosome names for consistency
    allelecounts$chr <- paste0("chr", allelecounts$chr)
  }
  
  # Create VRanges object
  allelecounts_vr <- tryCatch(
    {
      VRanges(
        seqnames = allelecounts$chr,
        ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos),
        ref = allelecounts$ref,
        alt = allelecounts$alt,
        seqinfo = seqinfo(bsgenome_selected)
      )
    },
    error = function(e) {
      stop(paste("Error creating VRanges object:", e$message))
    }
  )
  
  # Sort VRanges
  allelecounts_vr <- sort(allelecounts_vr)
  
  # Step 4: Assign Sample Name and Genotype Fields
  message("Assigning sample name and genotype fields...")

    # Assign sample name
  sampleNames(allelecounts_vr) <- sample_id
  
  # Create genotype data based on counts
  altDepth(allelecounts_vr) <- allelecounts$count_alt
  refDepth(allelecounts_vr) <- allelecounts$count_ref
  totalDepth(allelecounts_vr) <- allelecounts$count_ref + allelecounts$count_alt

  # **Add GT Field with '0/1' for Heterozygous Genotype**
  message("Adding 'GT' field with '0/1' for heterozygous genotype...")
  geno(allelecounts_vr)$GT <- matrix(
    rep("0/1", length(allelecounts_vr)),
    ncol = 1,
    dimnames = list(NULL, sample_id)
  )
  
  # Step 6: Write VCF to File
  message(paste("Writing VCF object to", outfile_vcf, "..."))
  tryCatch(
    {
      writeVcf(obj = allelecounts_vr, filename = outfile_vcf, index = TRUE)
      message("VCF file successfully written.")
    },
    error = function(e) {
      warning(paste("Error writing VCF file:", e$message))
    }
  )
  
  message("Function execution completed.")
  return(NULL)
}

countsdir <- "/staging/leuven/stg_00096/home/rdewin/ASE/results" 
sample_id <- "P011"
alias <- "tumor"


combine_loci_nomatch_enhanced(
  countsdir = countsdir,
  sample_id = sample_id,
  alias = alias
)



###TEST FUNCTION OF ASEREADCOUNT
ASEReadCount_improved <- function(sample_id,
                                   countsdir = "/staging/leuven/stg_00096/home/rdewin/ASE/results",
                                   RNA_dir = "/staging/leuven/stg_00096/home/rdewin/RNA/results/star",
                                   WGS_dir = "/staging/leuven/stg_00096/home/rdewin/WGS/resources",
                                   gatk_exe = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/gatk",
                                   java_cmd = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/bin/java",
                                   gatk_jar = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar",
                                   java_home = "/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/WGS",
                                   outfile_dir = NULL,
                                   minBaseQ = 20, minMapQ = 35) {
    # Load necessary libraries
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The 'readr' package is required but not installed.")
  }
  
  # Set JAVA_HOME environment variable
  Sys.setenv(JAVA_HOME = java_home)
  
  # Construct file paths based on sample_id and countsdir
  hetSNPvcf <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/P011_test.vcf.gz"
  #hetSNPvcf <- file.path(countsdir, sample_id, paste0(sample_id, "_hetSNPs_nomatch.vcf.bgz"))
  bamfile <- file.path(RNA_dir, sample_id, paste0(sample_id, "_Aligned.sortedByCoord.withRG.bam"))
  refgenome <- file.path(WGS_dir, "genome.fa")
  
  # Define output directory and file
  if (is.null(outfile_dir)) {
    outfile_dir <- file.path(countsdir, sample_id)
  }
  outfile <- file.path(outfile_dir, paste0(sample_id, "_asereadcounts_nomatch.tsv"))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outfile_dir)) {
    dir.create(outfile_dir, recursive = TRUE)
    message(paste("Created output directory:", outfile_dir))
  }
  
  # Check existence of input files and executables
  required_files <- list(
    "HetSNPVCF" = hetSNPvcf,
    "BAM file" = bamfile,
    "Reference Genome" = refgenome,
    "GATK JAR" = gatk_jar,
    "Java Executable" = java_cmd
  )
  
  missing_files <- names(required_files)[!file.exists(unlist(required_files))]
  
  if (length(missing_files) > 0) {
    stop(paste("The following required files/executables are missing:", 
               paste(missing_files, collapse = ", ")))
  }

  # Unzip the VCF file if it is compressed
  # if (grepl("\\.bgz$", hetSNPvcf)) {
  #   message("Unzipping the input VCF file...")
  #   system(paste("bgzip -d -c", shQuote(hetSNPvcf), ">", shQuote(sub(".bgz$", "", hetSNPvcf))))
  #   hetSNPvcf <- sub(".bgz$", "", hetSNPvcf)
  # }
  
  # Construct the ASEReadCounter command
  cmd <- paste0(
    shQuote(java_cmd), " -Xmx12G -jar ", shQuote(gatk_jar), " ASEReadCounter",
    " -I ", shQuote(bamfile),
    " -V ", shQuote(hetSNPvcf),
    " -R ", shQuote(refgenome),
    " -O ", shQuote(outfile),
    " --min-depth 1",
    " --min-mapping-quality ", minMapQ,
    " --min-base-quality ", minBaseQ,
    " --lenient",
    " 2> ", shQuote(file.path(outfile_dir, "ASEReadCounter.stderr"))
  )
  
  # Inform the user about the command being executed
  message("Executing ASEReadCounter command:")
  message(cmd)
  
  # Execute the command and capture the exit status
  exit_status <- system(cmd, wait = TRUE)
  
  # Check if the command was successful
  if (exit_status == 0) {
    message(paste("ASEReadCounter successfully executed. Output written to:", outfile))
  } else {
    warning("ASEReadCounter command failed. Please check the error messages above.")
  }
  
  # Return the path to the output file
  return(outfile)
}

ASEReadCount_improved(sample_id = "P011")
