import pandas as pd
import sys
import datetime

#logging
sys.stderr = open(snakemake.log[0], "w")
print(datetime.datetime.now(), file=sys.stderr)

sample_info = snakemake.params["sample_info"]

# List of strandedness for each sample
strandedness_list = sample_info['stranded'].tolist()

# Paths to the count files for each sample
paths = snakemake.input

def aggregate_counts(files, strandedness_list):
    all_counts = []
    for file, strand in zip(files, strandedness_list):
        df = pd.read_csv(file, sep="\t", header=None, usecols=[0, 2, 3], skiprows=4)
        df.columns = ["gene_id", "forward", "reverse"]
        relevant_column = "forward" if strand == "forward" else "reverse"
        count_data = df.set_index("gene_id")[relevant_column]
        count_data.name = file.split('/')[-2]  # Naming series after sample ID
        all_counts.append(count_data)

    combined_counts = pd.concat(all_counts, axis=1)
    combined_counts.to_csv(snakemake.output["counts"], sep="\t", header=True)

aggregate_counts(paths, strandedness_list)