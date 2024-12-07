.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hs1)

# Load the BSgenome
bsgenome_selected <- BSgenome.Hsapiens.UCSC.hs1
message("Chromosome names detected with 'chr' prefix. Using BSgenome.Hsapiens.UCSC.hs1.")


sample_id <- "P011"
# Read data
data <- read.table("/staging/leuven/stg_00096/home/rdewin/ASE/results/P011/P011_hetSNPs_nomatch.txt", header = TRUE, stringsAsFactors = FALSE)


# Create GRanges object
gr <- GRanges(seqnames = data$chr, ranges = IRanges(start = data$pos, width = 1))

# Set seqinfo of GRanges to match the genome
seqlevelsStyle(gr) <- seqlevelsStyle(bsgenome_selected)
seqinfo(gr) <- seqinfo(bsgenome_selected)[seqlevels(gr)]
# Create GRanges object


# Create REF and ALT
ref <- DNAStringSet(data$ref)
alt <- CharacterList(as.list(data$alt))

#alt <- as(alt, "SimpleList")

# Create QUAL and FILTER
qual <- rep(NA_real_, nrow(data))
filter <- rep("PASS", nrow(data))

# Create fixed DataFrame
fixed <- DataFrame(REF = ref, ALT = alt, QUAL = qual, FILTER = filter)

# Create genotype matrices
gt <- matrix("0/1", nrow = nrow(data), ncol = 1)
colnames(gt) <- sample_id

ad_array <- array(dim = c(nrow(data), 1, 2))
ad_array[,,1] <- as.integer(data$count_ref)
ad_array[,,2] <- as.integer(data$count_alt)
dimnames(ad_array) <- list(NULL, sample_id, c("Ref", "Alt"))

dp <- matrix(as.integer(data$count_ref + data$count_alt), nrow = nrow(data), ncol = 1)
colnames(dp) <- sample_id

# Create geno SimpleList
geno_list <- SimpleList(
  GT = gt,
  AD = ad_array,
  DP = dp
)

# Create colData
col_data <- DataFrame(row.names = sample_id)

# Create VCF object
vcf_obj <- VCF(
  rowRanges = gr,
  colData = col_data,
  fixed = fixed,
  geno = geno_list
)

# Create meta-information lines
# meta_info <- DataFrame(
#   row.names = c("fileformat", "phasing", "source"),
#   Value = c(
#     "VCFv4.1",
#     "unphased",
#     paste("VarkakaiantAnnotation", packageVersion("VariantAnnotation"))
#   )
# )

Create meta-information lines without 'fileDate' and 'fileformat'
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
writeVcf(vcf_obj, "output.vcf", index = FALSE)



meta_info <- DataFrame(
  row.names = c("fileformat", "phasing", "source", "fileDate"),
  Value = c(
    "VCFv4.3", 
    "unphased", 
    paste("VariantAnnotation", packageVersion("VariantAnnotation")),
    format(Sys.Date(), "%Y%m%d")  # fileDate comes last
  )
)

# Create FORMAT header (this part stays the same)
format_df <- DataFrame(
  Number = c("1", "G", "1"),
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

temp_vcf <- tempfile(fileext = ".vcf")
writeVcf(vcf_obj, temp_vcf, index = FALSE)
new_vcf <- readVcf(temp_vcf)

allelecounts_500 <- allelecounts[1:500,]
  # Step 3: Create VRanges Object
  message("Creating VRanges object...")
  
  # Determine if chromosome names start with "chr"
  sample_size <- min(100, nrow(allelecounts_500))
  chromosome_has_prefix <- any(
    grepl(
      pattern = "^chr",
      x = allelecounts_500$chr[sample(x = 1:nrow(allelecounts_500), size = sample_size, replace = TRUE)]
    )
  )
  
  # Select appropriate BSgenome based on chromosome naming
  if (chromosome_has_prefix) {
    bsgenome_selected <- BSgenome.Hsapiens.UCSC.hs1
    message("Chromosome names detected with 'chr' prefix. Using BSgenome.Hsapiens.UCSC.hs1.")
  } else {
    bsgenome_selected <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    message("Chromosome names detected without 'chr' prefix. Using BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0.")
    # Add "chr" prefix to chromosome names
    allelecounts_500$chr <- paste0("chr", allelecounts_500$chr)
  }
  
  # Create VRanges object
  allelecounts_500_vr <- tryCatch(
    {
      VRanges(
        seqnames = allelecounts_500$chr,
        ranges = IRanges(start = allelecounts_500$pos, end = allelecounts_500$pos),
        ref = allelecounts_500$ref,
        alt = allelecounts_500$alt,
        seqinfo = seqinfo(bsgenome_selected)
      )
    },
    error = function(e) {
      stop(paste("Error creating VRanges object:", e$message))
    }
  )
  
  # Sort VRanges
  allelecounts_500_vr <- sort(allelecounts_500_vr)
  
  # Assign sample name
  sampleNames(allelecounts_500_vr) <- sample_id


# Step 1: Prepare Genotype Fields (AD, DP, FT)

# AD: Combine ref and alt counts into a single column as "ref,alt"
AD_combined <- apply(cbind(allelecounts_500$count_ref, allelecounts_500$count_alt), 1, function(x) paste(x[1], x[2], sep = ","))
AD_matrix <- matrix(AD_combined, ncol = 1, dimnames = list(NULL, sample_id))

# DP: Total depth is the sum of ref and alt counts
DP_matrix <- matrix(allelecounts_500$count_ref + allelecounts_500$count_alt, ncol = 1, dimnames = list(NULL, sample_id))

# FT: Filter status, default is "."
FT_matrix <- matrix(rep(".", nrow(allelecounts_500)), ncol = 1, dimnames = list(NULL, sample_id))

# Step 2: Create the geno_list
geno_list <- SimpleList(
    AD = AD_matrix,  # Allelic depths for ref and alt alleles
    DP = DP_matrix,  # Total depth (sum of ref and alt)
    FT = FT_matrix   # Filter status (default ".")
)

# Step 3: Create the VCF object

# Make sure allelecounts_vr and col_data_vcf are already correctly created
vcf <- VCF(
    rowRanges = allelecounts_500_vr,    # Genomic ranges (chromosome, position)
    colData = col_data_vcf,         # Sample information (sample_id)
    geno = geno_list                # Genotype data (AD, DP, FT)
)

output_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/results"

# Step 4: Write the VCF file to disk
outfile_vcf <- file.path(output_dir, sample_id, paste0(sample_id, "test.vcf"))
writeVcf(vcf, outfile_vcf, index = TRUE)




.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hs1)

sample_id <- "P011"

# Read data
data <- read.table("/staging/leuven/stg_00096/home/rdewin/ASE/results/P011/P011_hetSNPs_nomatch.txt", header = TRUE, stringsAsFactors = FALSE)

# Load the BSgenome
bsgenome_selected <- BSgenome.Hsapiens.UCSC.hs1
message("Chromosome names detected with 'chr' prefix. Using BSgenome.Hsapiens.UCSC.hs1.")

# Create GRanges object
gr <- GRanges(seqnames = data$chr, ranges = IRanges(start = data$pos, width = 1))

# Set seqinfo of GRanges to match the genome
seqlevelsStyle(gr) <- seqlevelsStyle(bsgenome_selected)
seqinfo(gr) <- seqinfo(bsgenome_selected)[seqlevels(gr)]

# Create REF and ALT
ref <- DNAStringSet(data$ref)
alt <- CharacterList(as.list(data$alt))

# Create QUAL and FILTER
qual <- rep(NA_real_, nrow(data))
filter <- rep("PASS", nrow(data))

# Create fixed DataFrame
fixed <- DataFrame(REF = ref, ALT = alt, QUAL = qual, FILTER = filter)

# Create genotype matrices
gt <- matrix("0/1", nrow = nrow(data), ncol = 1)
colnames(gt) <- sample_id

ad_array <- array(dim = c(nrow(data), 1, 2))
ad_array[,,1] <- as.integer(data$count_ref)
ad_array[,,2] <- as.integer(data$count_alt)
dimnames(ad_array) <- list(NULL, sample_id, c("Ref", "Alt"))

dp <- matrix(as.integer(data$count_ref + data$count_alt), nrow = nrow(data), ncol = 1)
colnames(dp) <- sample_id

# Create geno SimpleList
geno_list <- SimpleList(
  GT = gt,
  AD = ad_array,
  DP = dp
)

# Create colData
col_data <- DataFrame(row.names = sample_id)

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
  Number = c("1", "R", "1"),
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

# Assign the header to the VCF object
header(vcf_obj) <- vcf_header

# Create VCF object
vcf_obj <- VCF(
  rowRanges = gr,
  colData = col_data,
  fixed = fixed,
  geno = geno_list
)

# Write VCF to a temporary file
temp_vcf <- tempfile()
writeVcf(vcf_obj, temp_vcf)

# Read the VCF file into R as lines
vcf_lines <- readLines(temp_vcf)

# Identify header lines and the rest of the content
header_indices <- grep("^##", vcf_lines)
column_header_index <- grep("^#CHROM", vcf_lines)

# Extract header lines
header_lines <- vcf_lines[header_indices]

# Remove duplicate 'fileDate' lines
header_lines <- header_lines[!grepl("^##fileDate=", header_lines)]

# Remove any existing 'fileformat' lines
header_lines <- header_lines[!grepl("^##fileformat=", header_lines)]

# Create the 'fileformat' line and place it at the beginning
fileformat_line <- "##fileformat=VCFv4.1"
header_lines <- c(fileformat_line, header_lines)

# Extract the rest of the VCF content (column header and data lines)
vcf_content <- c(header_lines, vcf_lines[column_header_index:length(vcf_lines)])

# Write the corrected VCF content to the final output file
writeLines(vcf_content, "output.vcf")
