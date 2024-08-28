#!/bin/bash

# Define variables
SAMPLE="P019"
CACHE_DIR="resources/vep/cache"
FASTA="resources/genome.fa"
FAI="resources/genome.fa.fai"
VCF_INPUT="results/mutect2/${SAMPLE}/${SAMPLE}_pass_variants.vcf"
VCF_OUTPUT="results/vep/${SAMPLE}/${SAMPLE}_mutect_annotated_all.vcf"
STATS_OUTPUT="results/vep/${SAMPLE}/${SAMPLE}_mutect_annotated_all.html"
GNOMAD_VCF="resources/gnomad.vcf.gz"
DBSNP_VCF="resources/dbSNPv155_common.vcf.gz"
CLINVAR_VCF="resources/clinvar.vcf.gz"
LOG="logs/vep/${SAMPLE}/mutect_annotated.log"

# Create directories if they do not exist
mkdir -p $(dirname ${VCF_OUTPUT})
mkdir -p $(dirname ${LOG})

# Run VEP annotation
echo "Running VEP annotation..."
vep \
    --input_file ${VCF_INPUT} \
    --output_file ${VCF_OUTPUT} \
    --dir_cache ${CACHE_DIR} \
    --fasta ${FASTA} \
    --format vcf \
    --vcf \
    --stats_file ${STATS_OUTPUT} \
    --species homo_sapiens_gca009914755v4 \
    --assembly T2T-CHM13v2.0 \
    --cache_version 106 \
    --custom file=${GNOMAD_VCF},short_name=gnomAD,format=vcf,fields=INFO \
    --custom file=${DBSNP_VCF},short_name=dbSNP,format=vcf,fields=RS \
    --custom file=${CLINVAR_VCF},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --offline \
    &> ${LOG}

echo "VEP annotation completed for sample ${SAMPLE}"
