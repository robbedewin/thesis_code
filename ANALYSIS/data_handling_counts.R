# ===================================================================
# R Notebook: Integrating Responsive and Induction Failure Cohorts
# ===================================================================


.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")

# Install necessary packages if not already installed
#install.packages(c("dplyr", "readxl", "biomaRt", "stringr", "tibble"))

# Load the libraries
library(dplyr)
library(readxl)
library(biomaRt)
library(stringr)
library(tibble)

# -------------------------------
# Load Clinical Data and Filter for Responsive Patients
# -------------------------------

# Define the path to your Excel file containing the supplementary data from Polonen et al.
excel_file <- "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/Supplementary_Data_Polonen.xlsx"

# Read Supplementary Table 1 (ST1_Clinical_Data) and Table 3 (ST3_Sample_Annotations)
st1_clinical_data <- read_excel(excel_file, sheet = "ST1_Clinical_Data")
st3_sample_annotations <- read_excel(excel_file, sheet = "ST3_Sample_Annotations")

# Filter for TARGET cohort and exclude Induction Failure
responsive_patients <- st1_clinical_data %>%
  filter(`In.TARGET.Cohort.(n=.264).RNASeq.and.WES` == "yes" & 
         Event.Type != "Induction Failure")

# Extract the patient IDs
responsive_patients_ids <- responsive_patients$USI


# Extract the clinical data (ETP.STATUS) for responsive patients
responsive_clinical_data <- st1_clinical_data %>%
  dplyr::filter(USI %in% responsive_patients_ids) %>%
  dplyr::mutate(condition = "responsive", ETP_status = ifelse(is.na(ETP.STATUS), "Unknown", ETP.STATUS)) %>%
  dplyr::select(sample_id = USI, condition, ETP_status)  


# Write the responsive_clinical_data to a file
write.table(responsive_clinical_data, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/responsive_clinical_data.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# -------------------------------
# Combine Responsive and Induction Failure Clinical Data
# -------------------------------

# Load the Induction Failure (IF) clinical data
IF_data <- read.table("/staging/leuven/stg_00096/home/rdewin/RNA/config/units_stranded_ETPstatus_RNA.tsv", header = TRUE, sep = "\t")

# Process the IF clinical data
IF_clinical_data <- IF_data %>%
    dplyr::select(sample_id = sample_name, ETP_status = ETPstatus) %>%
    dplyr::mutate(condition = "IF", ETP_status = ifelse(ETP_status == "notETP", "Non-ETP", ETP_status)) %>%
    dplyr::select(sample_id, condition, ETP_status)


# Combine responsive and IF clinical data
all_clinical_data <- rbind(responsive_clinical_data, IF_clinical_data)

# Write the all_clinical_data 
write.table(all_clinical_data, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/all_clinical_data.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



# -------------------------------
# Load and Filter Responsive RNA Counts Data, and Map Ensembl Gene IDs to Gene Symbols
# -------------------------------

# Read responsive counts data dowloaded from https://www.synapse.org/Synapse:syn54032669/files/ 
counts_data <- read.csv("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/TALL_X01_counts.tsv", sep = "\t", row.names = 1)

# Filter counts data to include only responsive samples
responsive_counts_data <- counts_data[, colnames(counts_data) %in% responsive_patients_ids]

# Write the filtered counts data
write.table(responsive_counts_data, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/responsive_counts_data.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


# initialize biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Clean Ensembl IDs by removing version numbers (e.g., ".14")
clean_ensembl_ids <- str_remove(rownames(responsive_counts_data), "\\..*$")

# Fetch gene symbols
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = clean_ensembl_ids,
  mart = mart
)

# Inspect the fetched gene annotations
head(gene_annotations)

# Merge gene symbols with counts data
responsive_counts_data$gene_id <- gene_annotations$hgnc_symbol[match(clean_ensembl_ids, gene_annotations$ensembl_gene_id)]

# Remove genes without a gene symbol
responsive_counts_mapped <- responsive_counts_data %>%
  filter(!is.na(gene_id) & gene_id != "")

# Aggregate counts by gene_id by summing numeric columns
responsive_counts_final <- responsive_counts_mapped %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}")) %>%
  ungroup()


# Verify the aggregation
print(paste("Number of unique genes after aggregation:", nrow(responsive_counts_final)))

# Check for any remaining duplicates
any_duplicates <- any(duplicated(responsive_counts_final$gene_id))
print(paste("Are there any duplicate gene symbols after aggregation?", any_duplicates))  # Should return FALSE


# Set gene_id as row names
responsive_counts_final <- responsive_counts_final %>%
  column_to_rownames(var = "gene_id")



# -------------------------------
# Load and Process InductionFailure Counts Data
# -------------------------------

# Read IF counts data from TSV file
IF_counts_data <- read.table("/staging/leuven/stg_00096/home/rdewin/RNA/results/counts/combined_counts.tsv", header = TRUE, sep = "\t", row.names = 1)

# The IF counts data contains a lot of LOC IDs, genes that are not well annotated and do not have gene symbols, so we separate them. 
# Separate LOC genes (starting with "LOC") and Non-LOC genes
IF_LOC_counts <- IF_counts_data %>%
  filter(str_detect(rownames(.), "^LOC"))

IF_nonLOC_counts <- IF_counts_data %>%
  filter(!str_detect(rownames(.), "^LOC"))

# Check the number of LOC and Non-LOC genes
print(paste("Number of LOC genes:", nrow(IF_LOC_counts)))      # 14701
print(paste("Number of Non-LOC genes:", nrow(IF_nonLOC_counts)))  # 30580


# -------------------------------
# Combine Responsive and Induction Failure Counts Data
# -------------------------------

# Identify common gene symbols between responsive and IF Non-LOC cohorts
print(paste("Number of responsive genes :", nrow(responsive_counts_final)))  # 20395
print(paste("Number of IF Non-LOC genes :", nrow(IF_nonLOC_counts)))  # 30580
common_genes <- intersect(rownames(responsive_counts_final), rownames(IF_nonLOC_counts))
print(paste("Number of common genes:", length(common_genes)))  # 18452

# Subset responsive counts to common genes
responsive_counts_common <- responsive_counts_final[common_genes, ]

# Subset IF counts to common genes
IF_counts_common <- IF_nonLOC_counts[common_genes, ]


# Combine responsive and IF counts data by binding columns
combined_counts <- cbind(responsive_counts_common, IF_counts_common)

# Verify the combined counts data
print(dim(combined_counts))  # Should be (common_genes, 264 + IF_samples)

# Save the combined counts data to a TSV file
write.table(combined_counts, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/combined_counts_common_genes.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Save LOC genes separately
write.table(IF_LOC_counts, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/IF_LOC_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

