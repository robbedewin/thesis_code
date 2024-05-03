library(ASCAT)

cwd <- getwd()
# Define the output directory
output_dir <- file.path(cwd, snakemake@output[[1]])
# Check if the output directory exists, create it if it does not
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Open the log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Define the input files
sample_id <- snakemake@wildcards[["sample"]]
tumor_bam <- file.path(cwd, snakemake@input[["tumor_bam"]])
normal_bam <- file.path(cwd, snakemake@input[["normal_bam"]])
allelecounter <- snakemake@config[["allelecounter"]]
alleles_prefix <- snakemake@config[["alleles_prefix"]]
loci_prefix <- snakemake@config[["loci_prefix"]]
GC_file <- snakemake@config[["GC_file"]]
nthreads <- snakemake@threads

# Update the paths for output files to be within the output directory
tumourLogR_file <- file.path(output_dir, paste0(sample_id, "_tumor_LogR.txt"))
tumourBAF_file <- file.path(output_dir, paste0(sample_id, "_tumor_BAF.txt"))
normalLogR_file <- file.path(output_dir, paste0(sample_id, "_normal_LogR.txt"))
normalBAF_file <- file.path(output_dir, paste0(sample_id, "_normal_BAF.txt"))

# Create the alleleFrequencies directory
alleleFrequencies_dir <- file.path(output_dir, "alleleFrequencies")
dir.create(alleleFrequencies_dir, recursive = TRUE, showWarnings = TRUE)

# Temporarily change the working directory to alleleFrequencies_dir
setwd(alleleFrequencies_dir)

# Run ascat.prepareHTS
ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = paste0(sample_id, "_tumor"),
  normalname = paste0(sample_id, "_normal"),
  allelecounter_exe = allelecounter,
  alleles.prefix = alleles_prefix,
  loci.prefix = loci_prefix,
  gender = "XX",
  genomeVersion = "hg38", # Adjust as necessary
  nthreads = nthreads,
  tumourLogR_file = tumourLogR_file,
  tumourBAF_file = tumourBAF_file,
  normalLogR_file = normalLogR_file,
  normalBAF_file = normalBAF_file,
  chrom_names = c(1:22, "X")
)

# Change the working directory back to output_dir
setwd(output_dir)

# Load the ASCAT objects
ascat.bc <- ascat.loadData(
  Tumor_LogR_file = tumourLogR_file, 
  Tumor_BAF_file = tumourBAF_file, 
  Germline_LogR_file = normalLogR_file, 
  Germline_BAF_file = normalBAF_file, 
  gender = 'XX', 
  genomeVersion = "hg38"
)

# Plot raw data before and after correction
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = GC_file)
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")

# Perform segmentation
ascat.bc <- ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "")

# Run ASCAT and calculate metrics
ascat.output <- ascat.runAscat(ascat.bc, gamma=1, write_segments = TRUE)
QC <- ascat.metrics(ascat.bc, ascat.output)

# Save the ASCAT objects in the working directory
ascat_objects_path <- file.path(getwd(), paste0(sample_id, "_ASCAT_objects.Rdata"))
save(ascat.bc, ascat.output, QC, file = ascat_objects_path)


# Close logging
sink()
sink(type = "message")

