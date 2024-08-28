#!/bin/bash

# Define variables
SAMPLE="P018"
CACHE_DIR="resources/vep/cache"
FASTA="resources/genome.fa"
FAI="resources/genome.fa.fai"
VCF_INPUT="results/mutect2/${SAMPLE}/${SAMPLE}_pass_variants.vcf"
VCF_OUTPUT="results/vep/${SAMPLE}/${SAMPLE}_mutect_annotated.vcf"
STATS_OUTPUT="results/vep/${SAMPLE}/${SAMPLE}_mutect_annotated.html"
LOG="logs/vep/${SAMPLE}/mutect_annotated.log"

# Create directories if they do not exist
mkdir -p $(dirname ${VCF_OUTPUT})
mkdir -p $(dirname ${LOG})


# Run VEP annotation
echo "Running VEP annotation..."
vep \
    --input_file ${VCF_INPUT} \
    --output_file ${VCF_OUTPUT} \
    --cache \
    --dir_cache ${CACHE_DIR} \
    --fasta ${FASTA} \
    --format vcf \
    --vcf \
    --offline \
    --species homo_sapiens_gca009914755v4 \
    --stats_file ${STATS_OUTPUT} \
    --assembly T2T-CHM13v2.0 \
    --cache_version 106 \
    --no_sift \
    --no_polyphen \
    &> ${LOG}

echo "VEP annotation completed for sample ${SAMPLE}"
