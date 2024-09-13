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
library(S4Vectors)
library(GenomicRanges)
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
  reference_alleles$count_ref <- ifelse(reference_alleles$ref == "A", sample_allele_counts$count_A,
                                        ifelse(reference_alleles$ref == "C", sample_allele_counts$count_C, 
                                               ifelse(reference_alleles$ref == "G", sample_allele_counts$count_G, sample_allele_counts$count_T)))
  reference_alleles$count_alt <- ifelse(reference_alleles$alt == "A", sample_allele_counts$count_A,
                                        ifelse(reference_alleles$alt == "C", sample_allele_counts$count_C, 
                                               ifelse(reference_alleles$alt == "G", sample_allele_counts$count_G, sample_allele_counts$count_T)))

  
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
mindepth <- 5  # The minimum depth required for filtering the alleles
sample_id <- "P011"  # Unique identifier for the sample being analyzed
alias <- "tumor"  # Optional alias (e.g., "tumor" or "normal")

# Run the function
get_allele_counts_filtered(
  chr = chr,
  reference_alleles_dir = reference_alleles_dir,
  sample_allele_counts_dir = sample_allele_counts_dir,
  mindepth = mindepth,
  sample_id = sample_id,
  alias = alias
)

chromosomes <- c(1:22, "X")
for (chr in chromosomes) {
  # Call the get_allele_counts_filtered function for each chromosome
  get_allele_counts_filtered(chr = chr, 
                             reference_alleles_dir = reference_alleles_dir,
                             sample_allele_counts_dir = sample_allele_counts_dir,
                             mindepth = mindepth,
                             sample_id = sample_id,
                             alias = alias)
  
  # Optionally print a message to track progress
  print(paste("Processed chromosome", chr))
}

countsdir <- "/staging/leuven/stg_00096/home/rdewin/ASE/results"  # Directory where the combined allele counts will be stored
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
  
  # Assign genotype fields using SimpleList from S4Vectors
  geno(allelecounts_vr) <- SimpleList(
    DP = depth_matrix,
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

# Function to execute GATK ASEReadCounter
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


ase_annotate <- function(asedf) {
  library(biomaRt)
  library(GenomicRanges)
  
  # Define chromosomal regions in the format "chr:start-end"
  chromregions <- paste(asedf$contig, asedf$position, asedf$position, sep = ":", collapse = ",")
  
  # Connect to Ensembl BioMart
  ensembl <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  # Retrieve gene annotations overlapping the ASE loci
  annot <- tryCatch(
    {
      getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
        filters = "chromosomal_region",
        mart = ensembl,
        values = chromregions
      )
    },
    error = function(e) {
      warning("Error fetching annotations from BioMart.")
      return(NULL)
    }
  )
  
  if (is.null(annot)) {
    asedf$gene <- NA
    return(asedf)
  }
  
  # Create GRanges objects for ASE loci and gene annotations
  asegr <- GRanges(
    seqnames = asedf$contig,
    ranges = IRanges(start = asedf$position, end = asedf$position)
  )
  
  annotgr <- GRanges(
    seqnames = annot$chromosome_name,
    ranges = IRanges(start = annot$start_position, end = annot$end_position),
    gene = annot$external_gene_name
  )
  
  # Find overlaps between ASE loci and gene annotations
  overlaps <- findOverlaps(query = asegr, subject = annotgr)
  
  # Initialize gene column
  asedf$gene <- NA
  
  # Assign gene names to overlapping ASE loci
  overlap_hits <- unique(queryHits(overlaps))
  for (hit in overlap_hits) {
    genes <- unique(annot$external_gene_name[subjectHits(overlaps)[queryHits(overlaps) == hit]])
    asedf$gene[hit] <- paste(genes, collapse = ",")
  }
  
  return(asedf)
}


compute_pvals_nomatch <- function(toutdir, tsample, filtercutoff = 0.01, exclude_bad_snps = FALSE) {
  library(readr)
  library(VariantAnnotation)
  library(BiocGenerics)
  
  # Define file paths
  asecountsfile <- file.path(toutdir, paste0(tsample, "_asereadcounts_nomatch.rtable"))
  genomecountsfile <- file.path(toutdir, paste0(tsample, "_hetSNPs_nomatch.txt"))
  
  # Read ASE counts
  asecounts <- tryCatch(
    {
      read_tsv(file = asecountsfile, col_names = TRUE, col_types = "ciccciiiiiiii")
    },
    error = function(e) {
      warning(paste("Error reading ASE counts file:", asecountsfile))
      return(NULL)
    }
  )
  
  if (is.null(asecounts)) return(NULL)
  
  # Read genome counts
  genomecounts <- tryCatch(
    {
      read_tsv(file = genomecountsfile, col_names = c("chr", "pos", "ref", "alt", "refCountGenome", "altCountGenome"), col_types = "ciccii")
    },
    error = function(e) {
      warning(paste("Error reading genome counts file:", genomecountsfile))
      return(NULL)
    }
  )
  
  if (is.null(genomecounts)) return(NULL)
  
  # Merge ASE counts with genome counts
  asecounts$contig <- sub(pattern = "^chr", replacement = "", x = asecounts$contig)
  asedf <- merge(x = asecounts, y = genomecounts, by.x = c("contig", "position"), by.y = c("chr", "pos"))
  
  # Select relevant columns
  asedf <- asedf[, c("contig", "position", "refAllele", "altAllele", "refCountGenome", "altCountGenome", "refCount", "altCount")]
  
  # Compute filter using beta-binomial confidence intervals
  asedf$filter <- apply(
    X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")],
    MARGIN = 1,
    FUN = function(x) {
      pconfint <- qbeta(c(filtercutoff / 2, 1 - filtercutoff / 2), shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1)
      min(pconfint)
    }
  )
  
  # Compute p-values
  asedf$pval <- apply(
    X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")],
    MARGIN = 1,
    FUN = function(x) {
      betabinom.test.ab(q = x["refCount"], size = x["refCount"] + x["altCount"],
                        shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided")
    }
  )
  
  # Initialize adjusted p-values
  asedf$padj <- 1
  
  # Identify test-worthy SNPs
  is_testworthy <- asedf$filter <= filtercutoff
  
  if (exclude_bad_snps) {
    # Read problematic loci
    probloci <- tryCatch(
      {
        read_tsv(file = "/staging/leuven/stg_00096/home/rdewin/ASE/ReferenceFiles/battenberg_problem_loci/probloci_270415.txt.gz", col_types = "ci")
      },
      error = function(e) {
        warning("Error reading problematic loci file. Proceeding without excluding bad SNPs.")
        return(NULL)
      }
    )
    
    if (!is.null(probloci)) {
      isbadsnp <- paste0(asedf$contig, "_", asedf$position) %in% paste0(probloci$Chr, "_", probloci$Pos)
      is_testworthy <- is_testworthy & !isbadsnp
    }
  }
  
  # Adjust p-values for test-worthy SNPs
  asedf$padj[is_testworthy] <- p.adjust(asedf$pval[is_testworthy], method = "fdr")
  
  return(asedf)
}


plot_ase_manhattan <- function(asedf) {
  library(ggplot2)
  library(dplyr)
  
  # Ensure contig is a factor with ordered levels
  asedf$contig <- factor(asedf$contig, levels = c(1:22, "X"))
  
  # Compute cumulative position for Manhattan plot
  asedf <- asedf %>%
    arrange(contig, position) %>%
    group_by(contig) %>%
    mutate(cumulative_pos = position + max(position) * (as.numeric(contig) - 1)) %>%
    ungroup()
  
  # Define significance threshold
  sig_threshold <- -log10(0.05)
  
  # Create the Manhattan plot
  p <- ggplot(data = asedf, aes(x = cumulative_pos, y = -log10(pval))) +
    geom_point(aes(color = contig %% 2 == 0), alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("skyblue", "navy")) +
    geom_hline(yintercept = sig_threshold, color = "red", linetype = "dashed") +
    labs(x = "Chromosome", y = "-log10(p-value)", title = "ASE Manhattan Plot") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # Add chromosome labels
  chr_means <- asedf %>%
    group_by(contig) %>%
    summarize(mean_pos = mean(cumulative_pos))
  
  p <- p + scale_x_continuous(
    breaks = chr_means$mean_pos,
    labels = chr_means$contig
  )
  
  return(p)
}


get_ase_aml_cell <- function(SAMPLEID) {
  library(parallel)
  library(readr)
  library(VariantAnnotation)
  library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
  library(BSgenome.Hsapiens.UCSC.hs1)
  library(biomaRt)
  library(ggplot2)
  
  # Define paths and parameters
  ALLELECOUNTER <- "/staging/leuven/stg_00096/software/alleleCount/4.0.0-GCCcore-6.4.0/bin/alleleCounter"
  ALLELESDIR <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/ReferenceFiles/"
  JAVA <- "/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-2.el8.x86_64/jre/bin/java"
  
  EXOMEDIR <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/data/exomes/cell_lines/mapped/"
  RNADIR <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/data/RNA/cell_lines/mapped/"
  RNAREFGENOME <- "/staging/leuven/stg_00096/home/rdewin/projects/ascat/chm13/ReferenceFiles/hg19/ucsc.hg19.fasta"
  
  TOUTDIR <- file.path("/staging/leuven/stg_00096/home/rdewin/ASE/results", SAMPLEID)
  
  minMapQ <- 35
  minBaseQ <- 20
  
  NCORES <- 1
  
  # Create output directory
  if (!dir.exists(TOUTDIR)) {
    dir.create(path = TOUTDIR, recursive = TRUE)
    message(paste("Created output directory:", TOUTDIR))
  }
  
  # Locate BAM files
  TBAMFILE <- list.files(path = EXOMEDIR, pattern = paste0(SAMPLEID, ".*bam$"), recursive = TRUE, full.names = TRUE)
  TRNABAMFILE <- list.files(path = RNADIR, pattern = paste0(SAMPLEID, ".*bam$"), recursive = TRUE, full.names = TRUE)
  
  if (length(TBAMFILE) == 0) {
    warning(paste("No tumor BAM file found for sample:", SAMPLEID))
    return(NULL)
  }
  if (length(TRNABAMFILE) == 0) {
    warning(paste("No RNA BAM file found for sample:", SAMPLEID))
    return(NULL)
  }
  
  # Step 1: Get allele counts for tumor exome
  chrs <- c(1:22, "X")
  
  if (!file.exists(file.path(TOUTDIR, paste0(SAMPLEID, "_hetSNPs_nomatch.vcf.bgz")))) {
    # Parallel processing across chromosomes
    mcmapply(
      FUN = get_allele_counts_filtered,
      chr = chrs,
      MoreArgs = list(
        reference_alleles_dir = ALLELESDIR,
        sample_allele_counts_dir = EXOMEDIR,
        mindepth = minBaseQ,  # Adjust as needed
        sample_id = SAMPLEID,
        alias = "tumor"
      ),
      mc.cores = NCORES
    )
    
    # Combine filtered allele counts into VCF
    combine_loci_nomatch_1(
      countsdir = "/staging/leuven/stg_00096/home/rdewin/ASE/results",
      sample_id = SAMPLEID,
      alias = "tumor"
    )
  }
  
  # Step 2: Run ASEReadCounter on RNA BAM
  ase_output_file <- file.path(TOUTDIR, paste0(SAMPLEID, "_asereadcounts_nomatch.rtable"))
  
  if (!file.exists(ase_output_file)) {
    ASEReadCount(
      hetSNPvcf = file.path(TOUTDIR, paste0(SAMPLEID, "_hetSNPs_nomatch.vcf.bgz")),
      bamfile = TRNABAMFILE,
      refgenome = RNAREFGENOME,
      outfile = ase_output_file,
      minBaseQ = minBaseQ,
      minMapQ = minMapQ
    )
  }
  
  # Step 3: Compute p-values
  asedf <- compute_pvals_nomatch(
    toutdir = TOUTDIR,
    tsample = SAMPLEID,
    exclude_bad_snps = FALSE
  )
  
  if (is.null(asedf)) {
    warning(paste("No ASE data to process for sample:", SAMPLEID))
    return(NULL)
  }
  
  # Step 4: Gene Annotation
  asedf <- ase_annotate(asedf = asedf[asedf$pval <= 0.01, ])
  
  # Order the data frame
  asedf$contig <- factor(x = asedf$contig, levels = c(1:22, "X"))
  asedf <- asedf[order(asedf$contig, asedf$position), ]
  
  # Write the annotated ASE results
  write_tsv(
    x = asedf,
    path = file.path(TOUTDIR, paste0(SAMPLEID, "_ase_out.txt"))
  )
  message(paste("ASE results written to:", file.path(TOUTDIR, paste0(SAMPLEID, "_ase_out.txt"))))
  
  # Step 5: Plotting
  p4 <- plot_ase_manhattan(asedf = asedf)
  
  # Save the Manhattan plot
  ggsave(
    filename = file.path(TOUTDIR, paste0(SAMPLEID, "_manhattan.png")),
    plot = p4,
    dpi = 300,
    width = 10,
    height = 3
  )
  message(paste("Manhattan plot saved to:", file.path(TOUTDIR, paste0(SAMPLEID, "_manhattan.png"))))
  
  return(NULL)
}
