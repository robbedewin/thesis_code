import pandas as pd

# Define the path to your samples.tsv file
config = {"samples": "/staging/leuven/stg_00096/home/rdewin/WGS/config/samples_10samples.tsv",
            "units": "/staging/leuven/stg_00096/home/rdewin/WGS/config/units_10samples.tsv"}

units = pd.read_csv(
    config["units"],
    sep="\t", 
    comment="#",
)

# Read the samples.tsv file into a DataFrame (if needed)
samples = pd.read_csv(
    config["samples"],
    sep="\t", 
    comment="#",
).set_index(["sample_name", "unit_name"], drop=False).sort_index()

# Generate a list of unique sample names from the units DataFrame
unique_samples = list(set(u.sample_name for u in units.itertuples()))
unique_datatypes = list(set(u.datatype for u in units.itertuples()))
unique_aliases = list(set(u.alias for u in units.itertuples()))


# Create a new 'identifier' column by combining 'sample_name', 'datatype', and 'alias'
units['identifier'] = units['sample_name'] + '_' + units['datatype'] + '_' + units['alias']

# Create separate 'identifier_rna' and 'identifier_atac' columns for RNA and ATAC data without the alias
units.loc[units['datatype'] == 'rna', 'identifier_rna'] = units['sample_name'] + '_rna'
units.loc[units['datatype'] == 'atac', 'identifier_atac'] = units['sample_name'] + '_atac'

# Set 'identifier', 'identifier_rna', and 'identifier_atac' as the index
units.set_index(['identifier', 'identifier_rna', 'identifier_atac'], inplace=True)


def get_fastq_rna(wildcards):
    identifier = f"{wildcards.sample}_rna"
    try:
        fastqs = units.loc[(slice(None), identifier), :]  # Use a slice to match any value in the 'identifier' level
        return {"r1": fastqs["fq1"].values[0].strip(), "r2": fastqs["fq2"].values[0].strip()}
    except KeyError:
        raise ValueError(f"Missing RNA FASTQ data for identifier: {identifier}")

def get_fastq_dna(wildcards):
    identifier = f"{wildcards.sample}_dna_{wildcards.alias}"
    try:
        fastqs = units.loc[(identifier, slice(None)), :]  # Use a slice to match any value in the 'identifier_rna' level
        return {"r1": fastqs["fq1"].values[0].strip(), "r2": fastqs["fq2"].values[0].strip()}
    except KeyError:
        raise ValueError(f"Missing DNA FASTQ data for identifier: {identifier}")

def get_fastq_atac(wildcards):
    identifier = f"{wildcards.sample}_atac"
    try:
        fastqs = units.loc[(slice(None), slice(None), identifier), :]  # Use slices to match any value in the other levels
        return {"r1": fastqs["fq1"].values[0].strip(), "r2": fastqs["fq2"].values[0].strip()}
    except KeyError:
        raise ValueError(f"Missing ATAC FASTQ data for identifier: {identifier}")

print(units)

class Wildcards:
    def __init__(self, sample):
        self.sample = sample
        
        
wildcards = Wildcards(sample='P013')
get_fastq_rna(wildcards)
print(get_fastq_rna(wildcards))