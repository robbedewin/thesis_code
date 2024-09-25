# Load necessary libraries
library(maftools)
library(GenomicRanges)
library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(BSgenome)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(rtracklayer)
library(stringr)
library(dplyr)
library(tidyr)
library(GenomeInfoDb)
library(circlize)
library(grid)


# Define the directory containing MAF files
maf_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/vcf2maf/"

# List all MAF files
maf_files <- list.files(path = maf_dir, pattern = "_mutect_pass_variants_annotated.maf$", recursive = TRUE, full.names = TRUE)

# Read all MAF files into a list
maf_list <- lapply(maf_files, read.maf)

# Merge all MAF objects into one
maf_combined <- merge_mafs(maf_list)


# Define the directory containing ASCAT segmentation files
ascat_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/"

# List all segmentation files
seg_files <- list.files(path = ascat_dir, pattern = "_tumor.segments.txt$", recursive = TRUE, full.names = TRUE)

# Read all segmentation files into a list
ascat_list <- lapply(seg_files, function(file) {
  seg_data <- fread(file)
  # Ensure consistent sample naming
  seg_data$Sample <- basename(dirname(file))
  return(seg_data)
})

# Combine all segmentation data into one data frame
ascat_combined <- rbindlist(ascat_list)


# Define the directory containing GRIDSS VCF files
gridss_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss/"

# List all high-confidence somatic VCF files
sv_files <- list.files(path = gridss_dir, pattern = "_high_confidence_somatic.vcf$", recursive = TRUE, full.names = TRUE)

# Read all VCF files into a list
sv_list <- lapply(sv_files, function(file) {
  # Read VCF file
  vcf <- readVcf(file, genome = "CHM13v2.0")
  # Convert to GRanges object
  gr <- breakpointRanges(vcf)
  # Add sample information
  gr$Sample <- basename(dirname(file))
  return(gr)
})

# Combine all SV data into one GRanges object
sv_combined <- do.call(c, sv_list)


# Read the gene annotation GTF file
gtf_file <- "/staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf"
gtf <- rtracklayer::import(gtf_file)

# Filter for 'transcript' entries
transcripts_gtf <- gtf[gtf$type == "transcript"]

# Create a data frame from the GRanges object
transcripts_df <- as.data.frame(transcripts_gtf)

# Group by gene_id and gene_name to get gene-level coordinates
genes_df <- transcripts_df %>%
  group_by(seqnames, gene_id, gene_name) %>%
  summarise(
    start = min(start),
    end = max(end),
    strand = unique(strand)
  ) %>%
  ungroup()

# Create GRanges object for genes
genes_gr <- GRanges(
  seqnames = genes_df$seqnames,
  ranges = IRanges(start = genes_df$start, end = genes_df$end),
  strand = genes_df$strand,
  gene_id = genes_df$gene_id,
  gene_name = genes_df$gene_name
)


# Get the top mutated genes from the MAF data
top_genes <- getGeneSummary(maf_combined)
top_genes <- top_genes[order(-top_genes$MutatedSamples),]
top_genes_list <- head(top_genes$Hugo_Symbol, n = 20)  # Adjust n as needed


# Example list of known T-ALL genes (you can replace or extend this list)
tall_genes_list <- c("NOTCH1", "FBXW7", "PTEN", "CDKN2A", "CDKN2B", "TAL1", "LMO2", "TLX1", "TLX3")

genes_of_interest <- unique(c(top_genes_list, tall_genes_list))


# Subset MAF data to genes of interest
maf_subset <- subsetMaf(maf = maf_combined, genes = genes_of_interest, includeSyn = FALSE)


# Convert ASCAT data to GRanges
cnv_gr <- GRanges(seqnames = ascat_combined$chr,
                  ranges = IRanges(start = ascat_combined$startpos, end = ascat_combined$endpos),
                  nMajor = ascat_combined$nMajor,
                  nMinor = ascat_combined$nMinor,
                  Sample = ascat_combined$Sample)


# Harmonize seqlevelsStyle chr1 vs 1
seqlevelsStyle(cnv_gr) <- seqlevelsStyle(genes_gr)


# Annotate CNVs with gene names
cnv_overlaps <- findOverlaps(cnv_gr, genes_gr)

# Check if overlaps are found
if (length(cnv_overlaps) > 0) {
  # Annotate CNVs with gene names
  cnv_annotated <- cnv_gr[queryHits(cnv_overlaps)]
  mcols(cnv_annotated)$gene_name <- mcols(genes_gr)$gene_name[subjectHits(cnv_overlaps)]
  
  # Keep only CNVs in genes of interest
  cnv_annotated <- cnv_annotated[mcols(cnv_annotated)$gene_name %in% genes_of_interest]
} else {
  warning("No overlaps found between CNVs and genes.")
}


# Annotate SVs with gene names
sv_overlaps <- findOverlaps(sv_combined, genes_gr)
sv_annotated <- sv_combined[queryHits(sv_overlaps)]
mcols(sv_annotated)$gene_name <- genes_gr$gene_name[subjectHits(sv_overlaps)]

# Keep only SVs in genes of interest
sv_annotated <- sv_annotated[mcols(sv_annotated)$gene_name %in% genes_of_interest]

# Extract mutation data as a data frame
mutations_df <- maf_subset@data[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Tumor_Sample_Barcode")]

# Convert mutations_df to data.frame
mutations_df <- as.data.frame(mutations_df)

# Update Tumor_Sample_Barcode in mutations_df
mutations_df$Tumor_Sample_Barcode <- gsub("_tumor$", "", mutations_df$Tumor_Sample_Barcode)


# Prepare mutation data frame using explicit dplyr namespace
mutations_df <- mutations_df %>%
  dplyr::mutate(Mutation = dplyr::case_when(
    Variant_Classification %in% c("Missense_Mutation") ~ "Missense",
    Variant_Classification %in% c("Nonsense_Mutation") ~ "Nonsense",
    Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frameshift",
    Variant_Classification %in% c("Splice_Site") ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Mutation)

# Create mutation matrix without list columns
mutation_matrix <- mutations_df %>%
  pivot_wider(
    names_from = Tumor_Sample_Barcode,
    values_from = Mutation,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%
  column_to_rownames("Hugo_Symbol")

# Define CNV status
cnv_status <- as.data.frame(mcols(cnv_annotated))
cnv_status <- cnv_status %>%
  dplyr::mutate(CNV = case_when(
    nMajor + nMinor == 0 ~ "Homozygous_Deletion",
    nMajor + nMinor == 1 ~ "Hemizygous_Deletion",
    nMajor + nMinor >= 5 ~ "Amplification",
    TRUE ~ "Neutral"
  )) %>%
  dplyr::select(gene_name, Sample, CNV)

# Create CNV matrix without list columns
cnv_matrix <- cnv_status %>%
  filter(CNV != "Neutral") %>%
  pivot_wider(
    names_from = Sample,
    values_from = CNV,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%
  column_to_rownames("gene_name")

# Prepare SV data frame
sv_status <- as.data.frame(mcols(sv_annotated))
sv_status <- sv_status %>%
  dplyr::mutate(SV = "Structural_Variant") %>%
  dplyr::select(gene_name, Sample, SV)


# Create SV matrix without list columns
sv_matrix <- sv_status %>%
  pivot_wider(
    names_from = Sample,
    values_from = SV,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%
  column_to_rownames("gene_name")


# Get all unique gene names and sample names across all matrices
all_genes <- unique(c(rownames(mutation_matrix), rownames(cnv_matrix), rownames(sv_matrix)))
all_samples <- unique(c(colnames(mutation_matrix), colnames(cnv_matrix), colnames(sv_matrix)))

# Initialize an empty alterations matrix
alterations_matrix <- matrix("", nrow = length(all_genes), ncol = length(all_samples))
rownames(alterations_matrix) <- all_genes
colnames(alterations_matrix) <- all_samples

# Function to fill alterations_matrix
fill_alterations_matrix <- function(alterations_matrix, alteration_matrix) {
  for (gene in rownames(alteration_matrix)) {
    for (sample in colnames(alteration_matrix)) {
      value <- alteration_matrix[gene, sample]
      if (!is.na(value) && value != "") {
        if (alterations_matrix[gene, sample] == "") {
          alterations_matrix[gene, sample] <- value
        } else {
          alterations_matrix[gene, sample] <- paste(alterations_matrix[gene, sample], value, sep = ";")
        }
      }
    }
  }
  return(alterations_matrix)
}

# Fill alterations_matrix with mutation data
alterations_matrix <- fill_alterations_matrix(alterations_matrix, mutation_matrix)

# Fill with CNV data
alterations_matrix <- fill_alterations_matrix(alterations_matrix, cnv_matrix)

# Fill with SV data
alterations_matrix <- fill_alterations_matrix(alterations_matrix, sv_matrix)

# Replace empty strings with NA
alterations_matrix[alterations_matrix == ""] <- NA

# Convert to a matrix of character strings
alterations_matrix <- apply(alterations_matrix, c(1, 2), as.character)



# Define alteration functions
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#FFFFFF", col = NA))
  },
  # Define functions for each mutation type
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#377EB8", col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#E41A1C", col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#4DAF4A", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#984EA3", col = NA))
  },
  Other = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#FF7F00", col = NA))
  },
  # CNV types
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w, h * 0.33, y = y + h * 0.33, gp = gpar(fill = "#FFD700", col = NA))
  },
  Hemizygous_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w, h * 0.33, y = y - h * 0.33, gp = gpar(fill = "#1E90FF", col = NA))
  },
  Homozygous_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w, h * 0.33, y = y - h * 0.33, gp = gpar(fill = "#00008B", col = NA))
  },
  # Structural Variants
  Structural_Variant = function(x, y, w, h) {
    grid.circle(x, y, r = min(unit(w, "npc"), unit(h, "npc")) * 0.5, gp = gpar(fill = "#A65628", col = NA))
  }
)

# Define colors for legend
col <- c(
  "Missense" = "#377EB8",
  "Nonsense" = "#E41A1C",
  "Frameshift" = "#4DAF4A",
  "Splice_Site" = "#984EA3",
  "Other" = "#FF7F00",
  "Amplification" = "#FFD700",
  "Hemizygous_Deletion" = "#1E90FF",
  "Homozygous_Deletion" = "#00008B",
  "Structural_Variant" = "#A65628"
)

# Save the plot to a PDF file
pdf("oncoplot_TALL_samples.pdf", width = 12, height = 8)

# Create oncoprint
oncoPrint(
  alterations_matrix,
  alter_fun = alter_fun,
  col = col,
  alter_fun_is_vectorized = FALSE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  column_title = "Mutational Landscape of T-ALL Samples",
  heatmap_legend_param = list(
    title = "Alterations",
    at = names(col),
    labels = names(col)
  )
)

# Close the PDF device
dev.off()
