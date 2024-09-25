.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R/library")

library(VariantAnnotation)
library(GenomicRanges)


# Read data
data <- read.table("/staging/leuven/stg_00096/home/rdewin/ASE/results/P011/P011_hetSNPs_nomatch.txt", header = TRUE, stringsAsFactors = FALSE)

# Create GRanges object
gr <- GRanges(seqnames = data$chr, ranges = IRanges(start = data$pos, width = 1))

# Create REF and ALT
ref <- DNAStringSet(data$ref)
alt <- CharacterList(data$alt)



# Create QUAL and FILTER
qual <- rep(NA_real_, nrow(data))
filter <- rep("PASS", nrow(data))

# Create fixed DataFrame
fixed <- DataFrame(REF = ref, ALT = alt, QUAL = qual, FILTER = filter)

# Create genotype matrices
gt <- matrix("0/1", nrow = nrow(data), ncol = 1)
colnames(gt) <- "Sample1"

ad_array <- array(dim = c(nrow(data), 1, 2))
ad_array[,,1] <- as.integer(data$count_ref)
ad_array[,,2] <- as.integer(data$count_alt)
dimnames(ad_array) <- list(NULL, "Sample1", c("Ref", "Alt"))

dp <- matrix(as.integer(data$count_ref + data$count_alt), nrow = nrow(data), ncol = 1)
colnames(dp) <- "Sample1"

# Create geno SimpleList
geno_list <- SimpleList(
  GT = gt,
  AD = ad_array,
  DP = dp
)

# Create colData
col_data <- DataFrame(row.names = "Sample1")

# Create VCF object
vcf_obj <- VCF(
  rowRanges = gr,
  colData = col_data,
  fixed = fixed,
  geno = geno_list
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
vcf_header <- VCFHeader(samples = "Sample1", header = DataFrameList(FORMAT = format_df))

# Set header in VCF object
header(vcf_obj) <- vcf_header

# Write VCF file
writeVcf(vcf_obj, "output.vcf")

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

