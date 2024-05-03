library(ASCAT)

# Open the log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# Define the input and output files using relative paths
sample_id <- snakemake@wildcards[["sample"]]
tumor_bam <- snakemake@input[["tumor_bam"]]
normal_bam <- snakemake@input[["normal_bam"]]
allelecounter <- snakemake@config[["allelecounter"]]
alleles_prefix <- snakemake@config[["alleles_prefix"]]
loci_prefix <- snakemake@config[["loci_prefix"]]
GC_file <- snakemake@config[["GC_file"]]
nthreads <- snakemake@threads

tumourLogR_file <- paste0("./results/ascat/", sample_id, "_tumor_LogR.txt")
tumourBAF_file <- paste0("./results/ascat/", sample_id, "_tumor_BAF.txt")
normalLogR_file <- paste0("./results/ascat/", sample_id, "_normal_LogR.txt")
normalBAF_file <- paste0("./results/ascat/", sample_id, "_normal_BAF.txt")

plot_prefix <- paste0("./results/ascat/plots/", sample_id, "_")


# ASCAT analysis steps
ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = paste0(sample_id, "_tumor"),
  normalname = paste0(sample_id, "_normal"),
  allelecounter_exe = allelecounter,
  alleles.prefix = alleles_prefix,
  loci.prefix = loci_prefix,
  gender = "XX",
  genomeVersion = "hg38", # This would only matter when gender=XY
  nthreads = nthreads,
  tumourLogR_file = tumourLogR_file,
  tumourBAF_file = tumourBAF_file,
  normalLogR_file = normalLogR_file,
  normalBAF_file = normalBAF_file,
  chrom_names = c(1:22, "X")
)

ascat.bc <- ascat.loadData(Tumor_LogR_file = tumourLogR_file, Tumor_BAF_file = tumourBAF_file, Germline_LogR_file = normalLogR_file, Germline_BAF_file = normalBAF_file, gender = 'XX', genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = paste0(plot_prefix, "Before_correction_"))
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = GC_file)
ascat.plotRawData(ascat.bc, img.prefix = paste0(plot_prefix, "After_correction_"))
ascat.bc <- ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = plot_prefix)
ascat.output <- ascat.runAscat(ascat.bc, gamma=1, write_segments = TRUE)
QC <- ascat.metrics(ascat.bc, ascat.output)
save(ascat.bc, ascat.output, QC, file = snakemake@output[["ascat_objects"]])

