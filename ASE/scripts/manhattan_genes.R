library(dplyr)
library(readr)

# Directory containing the files
manhattan_dir <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/manhattan_genes"

# List all .tsv files
files <- list.files(manhattan_dir, pattern = "_manhattan_annotated_genes\\.tsv$", full.names = TRUE)

# Function to read and standardize column types
read_and_normalize <- function(file) {
  read_tsv(file) %>%
    mutate(
      contig = as.character(contig)  # Ensure `contig` is character
    )
}

# Apply the function to all files and combine
combined_df <- bind_rows(lapply(files, read_and_normalize))

aggregated_df <- combined_df %>%
  group_by(gene) %>%
  summarise(
    num_samples = n_distinct(sample_id),
    contigs = paste(unique(contig), collapse = ", "),
    positions = paste(unique(position), collapse = ", "),
    samples = paste(unique(sample_id), collapse = ", "),    
    pval_min = min(pval, na.rm = TRUE),
    padj_min = min(padj, na.rm = TRUE)
  ) %>%
  arrange(gene)

# Write the result to a file
output_combined_file <- file.path(manhattan_dir, "combined_manhattan_annotated_genes.tsv")
write_tsv(aggregated_df, output_combined_file)

# Create and save the second file: ordered by sample count
aggregated_df_by_samples <- aggregated_df %>%
  arrange(desc(num_samples), gene)  # Arrange by number of samples and then by gene

output_combined_file_by_samples <- file.path(manhattan_dir, "combined_manhattan_annotated_genes_by_samples.tsv")
write_tsv(aggregated_df_by_samples, output_combined_file_by_samples)
