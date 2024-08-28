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
library(VGAM)

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


# Function to get alleles and filter based on minimum depth
get_alleles_chr_nomatch <- function(chr, allelesdir, countsdir, mindepth = 3, sample_id, alias = "tumor") {
  # Define the conversion map
  conversion <- c("1" = "A", "2" = "C", "3" = "G", "4" = "T")
  
  # Construct file paths
  allelesfile <- file.path(allelesdir, paste0("allele_T2T_chr", chr, ".txt"))
  countsfile <- file.path(countsdir, paste0(sample_id, "/alleleFrequencies/", sample_id, "_", alias, "_alleleFrequencies_chr", chr, ".txt"))
  outfile_allelecounts <- file.path(countsdir, paste0(sample_id, "_", alias, "_inform_alleles_nomatch_chr", chr, ".txt"))
  
  # Ensure the output directory exists
  outdir <- dirname(outfile_allelecounts)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Read the alleles data with original column names
  alleles <- read_tsv(file = allelesfile, col_names = TRUE, col_types = "iii")
  
  # Rename the columns to c("pos", "ref", "alt")
  colnames(alleles) <- c("pos", "ref", "alt")
  
  # Convert numeric alleles to nucleotide bases
  alleles$ref <- conversion[as.character(alleles$ref)]
  alleles$alt <- conversion[as.character(alleles$alt)]
  
  # Read the counts data with specified column types
  counts <- read_tsv(file = countsfile, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
  
  # Assign the chromosome name to the alleles data
  alleles <- cbind(chr = counts$chr, alleles)
  
  # Calculate reference and alternate allele counts
  alleles$count_ref <- ifelse(alleles$ref == "A", counts$count_A,
                              ifelse(alleles$ref == "C", counts$count_C, 
                                     ifelse(alleles$ref == "G", counts$count_G, counts$count_T)))
  alleles$count_alt <- ifelse(alleles$alt == "A", counts$count_A,
                              ifelse(alleles$alt == "C", counts$count_C, 
                                     ifelse(alleles$alt == "G", counts$count_G, counts$count_T)))
  
  # Filter data based on minimum depth
  output <- alleles[alleles$count_ref >= mindepth & alleles$count_alt >= mindepth, ]
  
  # Write the filtered data to output file
  write_tsv(x = output, path = outfile_allelecounts, col_names = F)
  
  return(NULL)
}

# This function combines allele counts and creates a VCF file using VRanges, handling chromosome name formats.
# This function combines allele counts and creates a VCF file using VRanges, handling chromosome name formats.
combine_loci_nomatch_1 <- function(countsdir, sample_id, alias="tumor") {
  # List files in the counts directory matching the pattern
  allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_nomatch_chr", full.names = TRUE)

  # Combine the allele counts files
  allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))
  
  # Write the combined allele counts to a file
  outfile <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  write_tsv(x = allelecounts, path = outfile, col_names = TRUE)
  
  # Write the combined allele counts to a VCF file
  locivcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))

  # Determine if the chromosome names start with "chr" and adjust if necessary
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

  # Create the geno matrices with the correct dimensions
  dp_matrix <- matrix(as.integer(allelecounts$count_ref + allelecounts$count_alt), nrow = length(allelecounts_vr), ncol = 1)
  ad_matrix <- matrix(c(allelecounts$count_ref, allelecounts$count_alt), nrow = length(allelecounts_vr), ncol = 2)

  # Ensure the column names of the matrices are set
  colnames(dp_matrix) <- sample_id
  colnames(ad_matrix) <- c(paste0(sample_id, "_R"), paste0(sample_id, "_A"))

  # Add genotype fields
  geno(allelecounts_vr) <- SimpleList(
    DP = dp_matrix,
    AD = ad_matrix
  )

  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  
  # Write the VRanges object to a VCF file
  writeVcf(obj = allelecounts_vr, filename = locivcf, index = TRUE)
  
  return(NULL)
}

# Execute the function
sample_id <- "P011"
countsdir <- "/staging/leuven/stg_00096/home/rdewin/ASE/"
combine_loci_nomatch_1(countsdir = countsdir, sample_id = sample_id)


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