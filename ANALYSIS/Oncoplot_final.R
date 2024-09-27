.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")
# Load necessary libraries
library(maftools)                    # For handling MAF files
library(GenomicRanges)               # For genomic range operations
library(data.table)                  # For efficient data manipulation
library(ComplexHeatmap)              # For creating complex heatmaps (oncoplot)
library(ggplot2)                     # For general plotting
library(BSgenome)                    # For genome data
library(VariantAnnotation)           # For reading VCF files
library(StructuralVariantAnnotation) # For SV annotation
library(rtracklayer)                 # For importing GTF files
library(stringr)                     # For string manipulation
library(dplyr)                       # For data manipulation
library(tidyr)                       # For data tidying
library(GenomeInfoDb)                # For genome information
library(circlize)                    # For circular visualization
library(grid)                        # For graphical objects
library(tibble)                      # For data frame operations


# ----------------------------- Explanation of script -----------------------------
# The script is designed to generate an oncoplot representing 
# the mutational landscape of T-cell Acute Lymphoblastic Leukemia 
# (T-ALL) samples by integrating mutation, copy number variation 
# (CNV), and structural variant (SV) data.

# Script Overview
# This script performs the following key steps:

# 1. Load Required Libraries: Imports all necessary R packages 
#    for data manipulation, genomic analyses, and visualization.

# 2. Read and Merge Mutation Data: Reads Mutation Annotation 
#    Format (MAF) files for each sample, merges them into a 
#    combined MAF object, and identifies the top mutated genes.

# 3. Read and Combine CNV Data: Reads ASCAT segmentation files 
#    for CNV analysis, combines them, and annotates CNVs with 
#    gene names.

# 4. Read and Combine SV Data: Reads GRIDSS VCF files for 
#    structural variants, combines them, and annotates SVs with 
#    gene names.

# 5. Prepare Gene Annotation Data: Imports gene annotations 
#    from a GTF file to map genomic coordinates to gene names.

# 6. Select Genes of Interest: Identifies genes of interest 
#    based on the top mutated genes and known T-ALL-related genes.

# 7. Subset and Process Mutation, CNV, and SV Data: Filters the 
#    mutation, CNV, and SV datasets to include only the genes 
#    of interest.

# 8. Create Matrices for Mutations, CNVs, and SVs: Converts the 
#    mutation, CNV, and SV data into matrices suitable for 
#    generating an oncoplot.

# 9. Combine Matrices into a Single Alterations Matrix: Integrates 
#    the mutation, CNV, and SV matrices into a single alterations 
#    matrix.

# 10. Generate the Oncoplot: Uses the ComplexHeatmap package to 
#     create an oncoplot that visually represents the mutational 
#     landscape across samples.

# ----------------------------- Step 1: Read Mutation Data -----------------------------

# Define the directory containing MAF files
maf_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/vcf2maf/"

# List all MAF files
maf_files <- list.files(
  path = maf_dir,
  pattern = "_mutect_pass_variants_annotated.maf$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all MAF files into a list
maf_list <- lapply(maf_files, read.maf)

# Merge all MAF objects into one
maf_combined <- merge_mafs(maf_list)

# ----------------------------- Step 2: Read CNV Data -----------------------------

# Define the directory containing ASCAT segmentation files
ascat_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/"

# List all segmentation files
seg_files <- list.files(
  path = ascat_dir,
  pattern = "_tumor.segments.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all segmentation files into a list and add sample names
ascat_list <- lapply(seg_files, function(file) {
  seg_data <- fread(file)
  # Ensure consistent sample naming by extracting from the directory name
  seg_data$Sample <- basename(dirname(file))
  return(seg_data)
})

# Combine all segmentation data into one data frame
ascat_combined <- rbindlist(ascat_list)

# ----------------------------- Step 3: Read SV Data -----------------------------

# Define the directory containing GRIDSS VCF files
gridss_dir <- "/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss/"

# List all high-confidence somatic VCF files
sv_files <- list.files(
  path = gridss_dir,
  pattern = "_high_confidence_somatic.vcf$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all VCF files into a list and add sample names
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

# ----------------------------- Step 4: Read Gene Annotation Data -----------------------------

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

# ----------------------------- Step 5: Identify Genes of Interest -----------------------------

# Get the top mutated genes from the MAF data
top_genes <- getGeneSummary(maf_combined)
top_genes <- top_genes[order(-top_genes$MutatedSamples), ]
top_genes_list <- head(top_genes$Hugo_Symbol, n = 20)  # Adjust n as needed

# List of known T-ALL genes (you can replace or extend this list)
tall_genes_list <- c("NOTCH1", "FBXW7", "PTEN", "CDKN2A", "CDKN2B", "TAL1", "LMO2", "TLX1", "TLX3")

# Combine top genes and known T-ALL genes
genes_of_interest <- unique(c(top_genes_list, tall_genes_list))

# ----------------------------- Step 6: Subset and Process Mutation Data -----------------------------

# Subset MAF data to genes of interest
maf_subset <- subsetMaf(maf = maf_combined, genes = genes_of_interest, includeSyn = FALSE)

# Extract mutation data as a data frame
mutations_df <- maf_subset@data[, c(
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Variant_Classification",
  "Tumor_Sample_Barcode"
)]

# Convert mutations_df to data.frame
mutations_df <- as.data.frame(mutations_df)

# Update Tumor_Sample_Barcode in mutations_df by removing '_tumor' suffix
mutations_df$Tumor_Sample_Barcode <- gsub("_tumor$", "", mutations_df$Tumor_Sample_Barcode)

# Map Variant_Classification to simplified Mutation types
mutations_df <- mutations_df %>%
  dplyr::mutate(Mutation = dplyr::case_when(
    Variant_Classification %in% c("Missense_Mutation") ~ "Missense",
    Variant_Classification %in% c("Nonsense_Mutation") ~ "Nonsense",
    Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frameshift",
    Variant_Classification %in% c("Splice_Site") ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Mutation)

# ----------------------------- Step 7: Subset and Process CNV Data -----------------------------

# Convert ASCAT data to GRanges
cnv_gr <- GRanges(
  seqnames = ascat_combined$chr,
  ranges = IRanges(start = ascat_combined$startpos, end = ascat_combined$endpos),
  nMajor = ascat_combined$nMajor,
  nMinor = ascat_combined$nMinor,
  Sample = ascat_combined$Sample
)

# Harmonize seqlevelsStyle (e.g., 'chr1' vs '1')
seqlevelsStyle(cnv_gr) <- seqlevelsStyle(genes_gr)

# Find overlaps between CNVs and genes
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

# Define CNV status based on copy number
cnv_status <- as.data.frame(mcols(cnv_annotated))
cnv_status <- cnv_status %>%
  dplyr::mutate(CNV = case_when(
    nMajor + nMinor == 0 ~ "Homozygous_Deletion",
    nMajor + nMinor == 1 ~ "Hemizygous_Deletion",
    nMajor + nMinor >= 5 ~ "Amplification",
    TRUE ~ "Neutral"
  )) %>%
  dplyr::select(gene_name, Sample, CNV)

# Filter out neutral CNVs
cnv_status <- cnv_status %>%
  filter(CNV != "Neutral")

# ----------------------------- Step 8: Subset and Process SV Data -----------------------------

# Find overlaps between SVs and genes
sv_overlaps <- findOverlaps(sv_combined, genes_gr)

# Annotate SVs with gene names
sv_annotated <- sv_combined[queryHits(sv_overlaps)]
mcols(sv_annotated)$gene_name <- genes_gr$gene_name[subjectHits(sv_overlaps)]

# Keep only SVs in genes of interest
sv_annotated <- sv_annotated[mcols(sv_annotated)$gene_name %in% genes_of_interest]

# Prepare SV data frame
sv_status <- as.data.frame(mcols(sv_annotated))
sv_status <- sv_status %>%
  dplyr::mutate(SV = "Structural_Variant") %>%
  dplyr::select(gene_name, Sample, SV)

# ----------------------------- Step 9: Create Matrices for Mutations, CNVs, and SVs -----------------------------

# Create mutation matrix
mutation_matrix <- mutations_df %>%
  pivot_wider(
    names_from = Tumor_Sample_Barcode,
    values_from = Mutation,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%  # Replace NA with empty strings
  column_to_rownames("Hugo_Symbol")

# Create CNV matrix
cnv_matrix <- cnv_status %>%
  pivot_wider(
    names_from = Sample,
    values_from = CNV,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%  # Replace NA with empty strings
  column_to_rownames("gene_name")

# Create SV matrix
sv_matrix <- sv_status %>%
  pivot_wider(
    names_from = Sample,
    values_from = SV,
    values_fn = function(x) paste(unique(x), collapse = ";")
  ) %>%
  replace(is.na(.), "") %>%  # Replace NA with empty strings
  column_to_rownames("gene_name")

# ----------------------------- Step 10: Combine Matrices into Alterations Matrix -----------------------------

# Get all unique gene names and sample names across all matrices
all_genes <- unique(c(
  rownames(mutation_matrix),
  rownames(cnv_matrix),
  rownames(sv_matrix)
))
all_samples <- unique(c(
  colnames(mutation_matrix),
  colnames(cnv_matrix),
  colnames(sv_matrix)
))

# Initialize an empty alterations matrix
alterations_matrix <- matrix(
  "",
  nrow = length(all_genes),
  ncol = length(all_samples),
  dimnames = list(all_genes, all_samples)
)

# Function to fill alterations_matrix with data
fill_alterations_matrix <- function(alterations_matrix, alteration_matrix) {
  common_genes <- intersect(rownames(alteration_matrix), rownames(alterations_matrix))
  common_samples <- intersect(colnames(alteration_matrix), colnames(alterations_matrix))
  
  for (gene in common_genes) {
    for (sample in common_samples) {
      value <- alteration_matrix[gene, sample]
      if (!is.na(value) && value != "") {
        if (alterations_matrix[gene, sample] == "") {
          alterations_matrix[gene, sample] <- value
        } else {
          alterations_matrix[gene, sample] <- paste(
            alterations_matrix[gene, sample],
            value,
            sep = ";"
          )
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

# Replace empty strings with NA for clarity
alterations_matrix[alterations_matrix == ""] <- NA

# Ensure all entries are character strings
alterations_matrix <- apply(alterations_matrix, c(1, 2), as.character)

# Replace NA values with empty strings to avoid warnings in oncoPrint
alterations_matrix[is.na(alterations_matrix)] <- ""

# ----------------------------- Step 11: Generate the Oncoplot -----------------------------

# Define alteration functions for oncoPrint
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

# Define colors for the legend
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
outputFilePath <- "/staging/leuven/stg_00096/home/rdewin/visualisation/oncoplot_TALL_samples_test2.pdf"
pdf(outputFilePath, width = 14, height = 14)

# Create the oncoprint
oncoPrint(
  alterations_matrix,
  alter_fun = alter_fun,
  col = col,
  alter_fun_is_vectorized = FALSE,  # Ensure correct rendering of multiple alterations
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


# Adjusted alter_fun with improved rect sizes and positions
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), 
              gp = gpar(fill = "#FFFFFF", col = "grey"))  # Add a border for clarity
  },
  # Mutation types
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w * 0.85, h * 0.85, gp = gpar(fill = "#377EB8", col = NA))  # Reduce size slightly
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w * 0.85, h * 0.85, gp = gpar(fill = "#E41A1C", col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w * 0.85, h * 0.85, gp = gpar(fill = "#4DAF4A", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w * 0.85, h * 0.85, gp = gpar(fill = "#984EA3", col = NA))
  },
  Other = function(x, y, w, h) {
    grid.rect(x, y, w * 0.85, h * 0.85, gp = gpar(fill = "#FF7F00", col = NA))
  },
  
  # CNV types with adjusted vertical positioning and reduced height
  Amplification = function(x, y, w, h) {
    grid.rect(x, y + h * 0.2, w, h * 0.2, gp = gpar(fill = "#FFD700", col = NA))  # Reduced height, adjusted y position
  },
  Hemizygous_Deletion = function(x, y, w, h) {
    grid.rect(x, y - h * 0.2, w, h * 0.2, gp = gpar(fill = "#1E90FF", col = NA))  # Reduced height, adjusted y position
  },
  Homozygous_Deletion = function(x, y, w, h) {
    grid.rect(x, y - h * 0.2, w, h * 0.2, gp = gpar(fill = "#00008B", col = NA))  # Reduced height, adjusted y position
  },
  
  # Structural Variants with circle shape to avoid overflow
  Structural_Variant = function(x, y, w, h) {
    grid.circle(x, y, r = min(unit(w, "npc"), unit(h, "npc")) * 0.4, gp = gpar(fill = "#A65628", col = NA))  # Circle with reduced radius
  }
)
