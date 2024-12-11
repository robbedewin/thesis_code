import pandas as pd

# Load VCF file
vcf_file = "/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss/P037/P037_high_confidence_somatic.vcf"
bnd_data = []

with open(vcf_file, "r") as file:
    for line in file:
        if not line.startswith("#"):
            fields = line.split("\t")
            chrom = fields[0]
            pos = fields[1]
            info = fields[7]
            mate = [x.split("=")[1] for x in info.split(";") if x.startswith("MATEID")]
            if mate:  # Check if MATEID is present
                bnd_data.append({"CHROM": chrom, "POS": pos, "MATEID": mate[0]})

bnd_df = pd.DataFrame(bnd_data)

# Debugging: Print the first few rows to check the MATEID field
print(bnd_df.head())

# Classify BNDs
def classify_bnd(row):
    mate_chrom = row["MATEID"].split(":")[0] if ":" in row["MATEID"] else None
    if mate_chrom:
        return "Inter-chromosomal" if row["CHROM"] != mate_chrom else "Intra-chromosomal"
    else:
        return "Unknown"

bnd_df["TYPE"] = bnd_df.apply(classify_bnd, axis=1)

# Summarize
summary = bnd_df["TYPE"].value_counts()
print(summary)