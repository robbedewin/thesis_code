sample_ids <- snakemake@params[['sample_ids']]

print(sample_ids)

# Write x to a text file
writeLines(sample_ids, con = "results/test/output.txt")
