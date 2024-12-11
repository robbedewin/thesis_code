# Load necessary libraries
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
library(ggplot2)


# Define the function to classify SVs
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "Inter-chromosomal", # inter-chromosomal
         ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "Insertion",
         ifelse(strand(gr) == strand(partner(gr)), "Inversion",
         ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "Deletion", "Duplication")))))
}

# Define the list of sample IDs
sample_ids <- c("P011", "P013", "P016", "P017", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P027", "P028", "P029", "P033", "P034", "P035", "P036", "P037", "P038", "P041", "P056", "P057", "P058", "P059", "P060", "P061", "P062", "P064", "P065", "P066", "P086", "P103", "P105")

# Initialize an empty dataframe to store SV type counts
sv_counts <- data.frame(Sample = character(), SVType = character(), Count = integer(), stringsAsFactors = FALSE)

# Loop through each sample and count SV types
for (sample_id in sample_ids) {
  # Load VCF file
  vcf_file <- file.path("/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss", sample_id, paste0(sample_id, "_high_confidence_somatic.vcf"))
  vcf <- readVcf(vcf_file, "CHM13_T2T")
  
  # Extract breakend (BND) records
  gr <- breakpointRanges(vcf)
  svtype <- simpleEventType(gr)
  
  # Count SV types
  svtype_counts <- table(svtype)
  
  # Add counts to the dataframe
  for (sv in names(svtype_counts)) {
    sv_counts <- rbind(sv_counts, data.frame(Sample = sample_id, SVType = sv, Count = svtype_counts[sv]))
  }
}

# Save the SV counts to a TSV file
write.table(sv_counts, "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/gridss/sv_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Summarize counts by SV type
sv_summary <- aggregate(Count ~ SVType, data = sv_counts, sum)

# Plot the distribution of SV types
ggplot(sv_summary, aes(x = SVType, y = Count, fill = SVType)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of Structural Variant Types", x = "SV Type", y = "Count") +
  scale_x_discrete(labels = c("Deletion", "Duplication", "Insertion", "Inversion", "Inter-chromosomal")) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/gridss/sv_type_distribution.png")

# Calculate averages and standard deviations per sample
sv_stats <- aggregate(Count ~ SVType, data = sv_counts, FUN = function(x) c(mean = mean(x), sd = sd(x)))
sv_stats <- do.call(data.frame, sv_stats)
colnames(sv_stats) <- c("SVType", "Mean", "SD")

# Plot the averages and standard deviations
ggplot(sv_stats, aes(x = SVType, y = Mean, fill = SVType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Average Number of Structural Variants per Sample", x = "SV Type", y = "Average Count") +
  scale_x_discrete(labels = c("Deletion", "Duplication", "Insertion", "Inversion", "Inter-chromosomal")) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/gridss/sv_type_averages.png")