#!/bin/bash

# Activate the Conda environment
source /lustre1/project/stg_00096/home/rdewin/system/miniconda/etc/profile.d/conda.sh
conda activate ASE_R

# Set the R_HOME environment variable
export R_HOME=/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/ASE_R/lib/R

# Path to your R script
R_SCRIPT_PATH="/staging/leuven/stg_00096/home/rdewin/ASE/scripts/ASE_final.R"

# Run the R script using Rscript
Rscript "$R_SCRIPT_PATH"