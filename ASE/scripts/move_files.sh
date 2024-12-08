#!/bin/bash

# Define the list of sample IDs
sample_ids=("P016" "P017" "P018" "P019" "P020" "P022" "P023" "P024" "P026" "P028" "P029" "P033" "P037" "P041" "P057" "P058" "P059" "P060" "P061" "P064" "P065" "P066" "P086" "P103" "P105")


# Define the base directories
base_dir="/staging/leuven/stg_00096/home/rdewin/ASE/results"
old_dir="${base_dir}/old"

# Loop through each sample ID and perform the operations
for sample_id in "${sample_ids[@]}"; do
  # Define the source and destination paths
  src_file1="${old_dir}/${sample_id}/${sample_id}_asereadcounts_nomatch.tsv"
  dest_file1="${base_dir}/${sample_id}/${sample_id}_asereadcounts.tsv"
  
  src_file2="${old_dir}/${sample_id}/ASEReadCounter.stderr"
  dest_file2="${base_dir}/${sample_id}/ASEReadCounter.stderr"
  
  src_folder="${old_dir}/${sample_id}/allele_counts"
  dest_folder="${base_dir}/${sample_id}"
  
  src_files3="${old_dir}/${sample_id}/${sample_id}_hetSNPs_nomatch.*"
  
  # Create the destination directory if it doesn't exist
  mkdir -p "${base_dir}/${sample_id}"
  
  # Copy the first file
  if [ -f "${src_file1}" ]; then
    cp "${src_file1}" "${dest_file1}"
    echo "Copied ${src_file1} to ${dest_file1}"
  else
    echo "File ${src_file1} does not exist"
  fi
  
  # Copy the second file
  if [ -f "${src_file2}" ]; then
    cp "${src_file2}" "${dest_file2}"
    echo "Copied ${src_file2} to ${dest_file2}"
  else
    echo "File ${src_file2} does not exist"
  fi
  
  # Move the folder
  if [ -d "${src_folder}" ]; then
    mv "${src_folder}" "${dest_folder}"
    echo "Moved ${src_folder} to ${dest_folder}"
  else
    echo "Folder ${src_folder} does not exist"
  fi
  
  # Copy the third set of files and remove the 'nomatch' label
  for file in ${src_files3}; do
    if [ -f "${file}" ]; then
      dest_file="${base_dir}/${sample_id}/$(basename "${file}" | sed 's/_nomatch//')"
      cp "${file}" "${dest_file}"
      echo "Copied ${file} to ${dest_file}"
    else
      echo "File ${file} does not exist"
    fi
  done
done