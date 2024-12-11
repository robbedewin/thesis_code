#!/bin/bash

# Initialize counters for each type of file
total_variants_type1=0
total_files_type1=0
total_variants_type2=0
total_files_type2=0

# Output file
output_file="ANALYSIS/mutect2/variant_counts.txt"
echo "Variant Counts Mutect" > "$output_file"


# Function to count variants in a VCF file
count_variants() {
    file=$1
    count=$(grep -v "^#" "$file" | wc -l)
    echo "$file: $count variants"
    if [[ "$file" == *"_filtered_variants.vcf" ]]; then
        total_variants_type1=$((total_variants_type1 + count))
        total_files_type1=$((total_files_type1 + 1))
    elif [[ "$file" == *"_pass_variants_new.vcf" ]]; then
        total_variants_type2=$((total_variants_type2 + count))
        total_files_type2=$((total_files_type2 + 1))
    fi
}

# Loop through the specified VCF files and count variants
for file in WGS/results/mutect2/*/*_filtered_variants.vcf WGS/results/mutect2/*/*_pass_variants_new.vcf; do
    if [[ -f "$file" ]]; then
        count_variants "$file"
    fi
done

# Calculate and print the average number of variants for each type
if [[ $total_files_type1 -gt 0 ]]; then
    average_variants_type1=$((total_variants_type1 / total_files_type1))
    echo "Average number of variants in *_filtered_variants: $average_variants_type1"
else
    echo "No *_filtered_variants files found."
fi

if [[ $total_files_type2 -gt 0 ]]; then
    average_variants_type2=$((total_variants_type2 / total_files_type2))
    echo "Average number of variants in *_pass_variants_new.vcf: $average_variants_type2"
else
    echo "No *_pass_variants_new.vcf files found."
fi