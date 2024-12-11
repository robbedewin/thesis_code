#!/bin/bash

# Define the list of sample IDs
sample_ids=("P011" "P013" "P016" "P017" "P018" "P019" "P020" "P022" "P023" "P024" "P026" "P027" "P028" "P029" "P033" "P034" "P035" "P036" "P037" "P038" "P041" "P056" "P057" "P058" "P059" "P060" "P061" "P062" "P064" "P065" "P066" "P086" "P103" "P105")

# Define the base directory
base_dir="/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss"

# Initialize counters and sums
count_all=0
sum_all=0
count_filtered=0
sum_filtered=0
count_high=0
sum_high=0

# Output file
output_file="ANALYSIS/gridss/variant_counts.txt"
echo "Variant Counts" > "$output_file"

# Loop through each sample ID and count the variants
for sample_id in "${sample_ids[@]}"; do
  # Define the file paths
  file_all="${base_dir}/${sample_id}/${sample_id}_all_calls.vcf"
  file_filtered="${base_dir}/${sample_id}/${sample_id}_high_and_low_confidence_somatic.vcf"
  file_high="${base_dir}/${sample_id}/${sample_id}_high_confidence_somatic.vcf"
  
  # Count the number of variants in each file
  if [ -f "${file_all}" ]; then
    num_all=$(grep -vc "^#" "${file_all}")
    sum_all=$((sum_all + num_all))
    count_all=$((count_all + 1))
    echo "${file_all}: ${num_all} variants" >> "$output_file"
  fi
  
  if [ -f "${file_filtered}" ]; then
    num_filtered=$(grep -vc "^#" "${file_filtered}")
    sum_filtered=$((sum_filtered + num_filtered))
    count_filtered=$((count_filtered + 1))
    echo "${file_filtered}: ${num_filtered} variants" >> "$output_file"
  fi
  
  if [ -f "${file_high}" ]; then
    num_high=$(grep -vc "^#" "${file_high}")
    sum_high=$((sum_high + num_high))
    count_high=$((count_high + 1))
    echo "${file_high}: ${num_high} variants" >> "$output_file"
  fi
done

# Calculate the averages
avg_all=$(echo "scale=2; $sum_all / $count_all" | bc)
avg_filtered=$(echo "scale=2; $sum_filtered / $count_filtered" | bc)
avg_high=$(echo "scale=2; $sum_high / $count_high" | bc)

# Print the averages at the top of the output file
{
  echo "Average number of variants (all): $avg_all"
  echo "Average number of variants (filtered high and low confidence): $avg_filtered"
  echo "Average number of variants (only high confidence): $avg_high"
  echo ""
  cat "$output_file"
} > temp_output_file && mv temp_output_file "$output_file"