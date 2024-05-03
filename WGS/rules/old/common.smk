import pandas as pd
import gzip
import os

# Read the units.tsv file into a DataFrame
units_path = "config/units_test.tsv"
units = pd.read_csv(units_path, sep="\t", dtype=str, header=0)

# Strip spaces from column names
units.columns = units.columns.str.strip()

# Print the unique sample and unit names
#print(units['sample_name'].unique())
#print(units['unit_name'].unique())

# Now set the index
units.set_index([units["sample_name"].str.strip(), units["unit_name"].str.strip()], inplace=True)

def get_fastq(wildcards):
    """Get fastq files for a given sample and unit."""
    # Access the row for the given sample and unit using the MultiIndex
    try:
        fastqs = units.loc[(wildcards.sample, wildcards.unit)]
        # Strip spaces from the fastq file paths and return them
        return {"r1": fastqs["fq1"].strip(), "r2": fastqs["fq2"].strip()}
    except KeyError:
        raise ValueError(f"Missing data for sample {wildcards.sample} and unit {wildcards.unit}")


def extract_read_group(fastq_path):
    # Check if the file is compressed and open accordingly
    open_func = gzip.open if fastq_path.endswith('.gz') else open
    with open_func(fastq_path, 'rt') as fastq_file:  # 'rt' mode for reading text
        first_line = fastq_file.readline().strip()
    parts = first_line.split(':')
    flowcell = parts[2]
    lane = parts[3]
    barcode = parts[-1]
    # Construct the read group string
    rg = f"@RG\\tID:{flowcell}.{lane}\\tPL:ILLUMINA\\tPU:{flowcell}.{lane}.{barcode}\\tSM:sample"
    return rg