# Load necessary libraries
library(circlize)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(data.table)
library(VariantAnnotation)
library(rtracklayer)

# Load MAF files
maf_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/vcf2maf/"
maf_files <- list.files(
  path = maf_dir,
  pattern = "_mutect_pass_variants_annotated.maf$",
  recursive = TRUE,
  full.names = TRUE
)

maf_list <- lapply(maf_files, function(file) {
  fread(file)
})

mutations_df <- do.call(rbind, maf_list)

# Load CNV files
ascat_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/"
seg_files <- list.files(
  path = ascat_dir,
  pattern = "_tumor.segments.txt$",
  recursive = TRUE,
  full.names = TRUE
)

cnv_list <- lapply(seg_files, function(file) {
  data <- fread(file)
  data$Sample <- basename(dirname(file))
  return(data)
})

cnv_status <- do.call(rbind, cnv_list)

# Load SV files
gridss_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss/"
sv_files <- list.files(
  path = gridss_dir,
  pattern = "_high_confidence_somatic.vcf$",
  recursive = TRUE,
  full.names = TRUE
)

sv_list <- lapply(sv_files, function(file) {
  vcf <- readVcf(file, genome = "CHM13v2.0")
  gr <- breakpointRanges(vcf)
  gr$Sample <- basename(dirname(file))
  return(gr)
})

sv_status <- do.call(c, sv_list)

# Load gene annotation data
gtf_file <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf"
genes_gtf <- rtracklayer::import(gtf_file)
genes_gr <- genes_gtf[genes_gtf$type == "gene"]

# Annotate CNV data
annotate_cnvs <- function(cnv_status, genes_gr) {
  # Convert CNV data to GRanges
  cnv_gr <- GRanges(
    seqnames = cnv_status$chr,
    ranges = IRanges(start = cnv_status$startpos, end = cnv_status$endpos),
    Sample = cnv_status$Sample
  )

  # Harmonize seqlevelsStyle
  seqlevelsStyle(cnv_gr) <- seqlevelsStyle(genes_gr)

  # Find overlaps between CNVs and gene annotations
  overlaps <- findOverlaps(cnv_gr, genes_gr)

  # Annotate CNVs with gene names
  cnv_annotated <- cnv_gr[queryHits(overlaps)]
  mcols(cnv_annotated)$gene_name <- mcols(genes_gr)$gene_name[subjectHits(overlaps)]

  # Convert to data frame with genomic coordinates
  cnv_df <- as.data.frame(mcols(cnv_annotated))
  cnv_df$seqnames <- as.character(seqnames(cnv_annotated))
  cnv_df$start <- start(cnv_annotated)
  cnv_df$end <- end(cnv_annotated)
  return(cnv_df)
}

cnv_status <- annotate_cnvs(cnv_status, genes_gr)

# Annotate SV data
annotate_svs <- function(sv_status, genes_gr) {
  # Find overlaps between SVs and gene annotations
  overlaps <- findOverlaps(sv_status, genes_gr)

  # Annotate SVs with gene names
  sv_annotated <- sv_status[queryHits(overlaps)]
  mcols(sv_annotated)$gene_name <- mcols(genes_gr)$gene_name[subjectHits(overlaps)]

  # Convert to data frame with genomic coordinates
  sv_df <- as.data.frame(mcols(sv_annotated))
  sv_df$seqnames <- as.character(seqnames(sv_annotated))
  sv_df$start <- start(sv_annotated)
  sv_df$end <- end(sv_annotated)
  return(sv_df)
}

sv_status <- annotate_svs(sv_status, genes_gr)

# Define function to prepare data for circos plot
prepare_circos_data <- function(mutations_df, cnv_status, sv_status) {
  # Combine mutation, CNV, and SV data
  mutations_circos <- mutations_df %>%
    mutate(type = "Mutation") %>%
    select(seqnames = Chromosome, start = Start_Position, end = End_Position, type)

  cnv_circos <- cnv_status %>%
    mutate(type = CNV) %>%
    select(seqnames, start, end, type)

  sv_circos <- sv_status %>%
    mutate(type = "SV") %>%
    select(seqnames, start, end, type)

  # Combine all alterations
  combined_data <- bind_rows(mutations_circos, cnv_circos, sv_circos)

  return(combined_data)
}

# Prepare circos data
circos_data <- prepare_circos_data(mutations_df, cnv_status, sv_status)

# Define colors for each alteration type
alteration_colors <- c(
  "Mutation" = "#377EB8",
  "Amplification" = "#FFD700",
  "Hemizygous_Deletion" = "#1E90FF",
  "Homozygous_Deletion" = "#00008B",
  "SV" = "#A65628"
)

# Initialize circos plot with genomic data
circos.initializeWithIdeogram(plotType = c("axis", "labels"))

# Add tracks for mutations, CNVs, and SVs
circos.trackPlotRegion(
  factors = circos_data$seqnames,
  y = seq_len(nrow(circos_data)),
  ylim = c(0, 1),
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    chr <- CELL_META$sector.index
    chr_data <- circos_data[circos_data$seqnames == chr, ]

    # Plot mutations
    mutation_data <- chr_data[chr_data$type == "Mutation", ]
    if (nrow(mutation_data) > 0) {
      circos.points(
        x = mutation_data$start,
        y = rep(0.5, nrow(mutation_data)),
        col = alteration_colors["Mutation"],
        pch = 16,
        cex = 0.6
      )
    }

    # Add CNVs
    cnv_data <- chr_data[chr_data$type %in% c("Amplification", "Hemizygous_Deletion", "Homozygous_Deletion"), ]
    if (nrow(cnv_data) > 0) {
      circos.rect(
        xleft = cnv_data$start,
        xright = cnv_data$end,
        ybottom = 0.1,
        ytop = 0.3,
        col = alteration_colors[cnv_data$type],
        border = NA
      )
    }

    # Add SVs
    sv_data <- chr_data[chr_data$type == "SV", ]
    if (nrow(sv_data) > 0) {
      circos.links(
        x1 = sv_data$start,
        x2 = sv_data$end,
        col = alteration_colors["SV"],
        border = NA,
        lwd = 1
      )
    }
  }
)

# Add legend
legend("topright", legend = names(alteration_colors), fill = alteration_colors, border = "black", bty = "n")

# Save plot to file
outputFilePath <- "/staging/leuven/stg_00096/home/rdewin/PLOTS/circos_plot.pdf"
pdf(outputFilePath, width = 8, height = 8)
circos.clear()
dev.off()
