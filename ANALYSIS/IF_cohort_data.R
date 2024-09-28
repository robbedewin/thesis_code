.libPaths("/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/R/lib/R/library")

library(dplyr)

# Load the data
IF_data <- read.table("/staging/leuven/stg_00096/home/rdewin/RNA/config/units_stranded_ETPstatus_RNA.tsv", header = TRUE, sep = "\t")



IF_clinical_data <- IF_data %>%
    select(sample_id = sample_name, ETP_status = ETPstatus) %>%
    mutate(condition = "IF", ETP_status = ifelse(ETP_status == "notETP", "Non-ETP", ETP_status)) %>%
    select(sample_id, condition, ETP_status)

# Write the IF_clinical_data
write.table(IF_clinical_data, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/IF_clinical_data.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# Combine responsive and IF clinical data
all_clinical_data <- rbind(responsive_clinical_data, IF_clinical_data)

# Write the all_clinical_data
write.table(all_clinical_data, file = "/staging/leuven/stg_00096/home/rdewin/ANALYSIS/data/all_clinical_data.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
