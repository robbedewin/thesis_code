# Define the directory containing the log files
log_dir <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/star"

# List all Log.final.out files
log_files <- list.files(path = log_dir, pattern = "_Log.final.out$", recursive = TRUE, full.names = TRUE)

# Define the matched samples
matchedSamples <- c("P011", "P013", "P016", "P017", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P028", "P029", "P033", "P037", "P041", "P057", "P058", "P059", "P060", "P061", "P064", "P065", "P066", "P086", "P103", "P105")

# Filter the log files to include only those with sample IDs in matchedSamples
filtered_log_files <- log_files[sapply(log_files, function(x) any(sapply(matchedSamples, function(y) grepl(y, x))))]

# Initialize vectors to store statistics
alignment_rates <- numeric()
unique_reads <- numeric()
multi_mapping_reads <- numeric()
unmapped_reads <- numeric()

# Function to extract statistics from a log file
extract_stats <- function(file) {
  lines <- readLines(file)
  alignment_rate <- as.numeric(sub(".*\\|\\s*(\\d+\\.\\d+)%$", "\\1", lines[grep("Uniquely mapped reads %", lines)]))
  unique_read_count <- as.numeric(sub(".*\\|\\s*(\\d+)$", "\\1", lines[grep("Uniquely mapped reads number", lines)]))
  multi_mapping_read_count <- as.numeric(sub(".*\\|\\s*(\\d+)$", "\\1", lines[grep("Number of reads mapped to multiple loci", lines)]))
  unmapped_read_count <- as.numeric(sub(".*\\|\\s*(\\d+)$", "\\1", lines[grep("Number of reads unmapped: too short", lines)]))
  
  return(list(alignment_rate = alignment_rate, unique_reads = unique_read_count, multi_mapping_reads = multi_mapping_read_count, unmapped_reads = unmapped_read_count))
}

# Loop through each filtered log file and extract statistics
for (file in filtered_log_files) {
  stats <- extract_stats(file)
  alignment_rates <- c(alignment_rates, stats$alignment_rate)
  unique_reads <- c(unique_reads, stats$unique_reads)
  multi_mapping_reads <- c(multi_mapping_reads, stats$multi_mapping_reads)
  unmapped_reads <- c(unmapped_reads, stats$unmapped_reads)
}

# Calculate averages
average_alignment_rate <- mean(alignment_rates)
standard_deviation_alignment_rate <- sd(alignment_rates)
average_unique_reads <- mean(unique_reads)
average_multi_mapping_reads <- mean(multi_mapping_reads)
average_unmapped_reads <- mean(unmapped_reads)

# Print the results
cat("Average Alignment Rate:", average_alignment_rate, "%\n")
cat("Standard Deviation Alignment Rate:", standard_deviation_alignment_rate, "\n")
cat("Average Unique Reads:", average_unique_reads, "\n")
cat("Average Multi-Mapping Reads:", average_multi_mapping_reads, "\n")
cat("Average Unmapped Reads:", average_unmapped_reads, "\n")
