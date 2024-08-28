import pandas as pd

# Define the path to your samples.tsv file
config = {"samples": "/staging/leuven/stg_00096/home/rdewin/WGS/config/samples_10samples.tsv",
            "units": "/staging/leuven/stg_00096/home/rdewin/RNA/config/units_10samples_stranded_ETPstatus.tsv"}

units = pd.read_csv(
    config["units"],
    sep="\t", 
    comment="#",
)

unique_samples = list(set(u.sample_name for u in units.itertuples()))

def get_strandedness(sample_name):
    subset = units.loc[(units["sample_name"] == sample_name) & (units["datatype"] == "rna"), "stranded"]
    return subset.values[0] if not subset.empty else "unknown"
    
def get_ETPstatus(sample_name):
    units = pd.read_csv(config["units"], sep="\t", comment="#")
    subset = units.loc[(units["sample_name"] == sample_name) & (units["datatype"] == "rna"), "ETPstatus"]
    return subset.values[0] if not subset.empty else "unknown"

etp_status_list = [get_ETPstatus(sample) for sample in units["sample_name"].unique()]

samples = units["sample_name"].unique()
sample_info = pd.DataFrame({
    "sample_name": samples,
    "stranded": [get_strandedness(sample) for sample in samples],
    "ETPstatus": [get_ETPstatus(sample) for sample in samples]
    
})

# Save the sample information to a CSV file
sample_info_file = "RNA/results/sample_info.csv"
sample_info.to_csv(sample_info_file, index=False)

etp_status_list = [get_ETPstatus(sample) for sample in units["sample_name"].unique()]



#print(get_strandedness('P027'))

for sample in unique_samples:
    print(sample, get_strandedness(sample), get_ETPstatus(sample))
    
 # Generate a list of unique sample names from the units DataFrame


# # Use the expand function
# paths = ["results/gridss/{}_normal.vcf".format(sample) for sample in unique_samples]

# # Print the results
# for path in paths:
#     print(path)