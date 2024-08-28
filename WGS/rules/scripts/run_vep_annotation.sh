#!/bin/bash

# This script runs VEP annotation with custom annotations
# Usage: ./run_vep_annotation.sh <sample> <input_vcf> <output_vcf> <output_html> <cache_dir> <fasta_file> <gnomad_vcf> <dbsnp_vcf> <clinvar_vcf> <log_file> <cache_version> <species> <assembly>

SAMPLE=$1
VCF_INPUT=$2
VCF_OUTPUT=$3
STATS_OUTPUT=$4
CACHE_DIR=$5
FASTA=$6
GNOMAD_VCF=$7
DBSNP_VCF=$8
CLINVAR_VCF=$9
LOG=${10}
CACHE_VERSION=${11}
SPECIES=${12}
ASSEMBLY=${13}

# Create directories if they do not exist
mkdir -p $(dirname ${VCF_OUTPUT})
mkdir -p $(dirname ${LOG})

echo "Running VEP annotation..."
vep \
    --input_file ${VCF_INPUT} \
    --output_file ${VCF_OUTPUT} \
    --offline \
    --cache \
    --dir_cache ${CACHE_DIR} \
    --fasta ${FASTA} \
    --format vcf \
    --vcf \
    --stats_file ${STATS_OUTPUT} \
    --species ${SPECIES} \
    --assembly ${ASSEMBLY} \
    --cache_version ${CACHE_VERSION} \
    --custom file=${GNOMAD_VCF},short_name=gnomAD,format=vcf,fields=INFO \
    --custom file=${DBSNP_VCF},short_name=dbSNP,format=vcf,fields=RS \
    --custom file=${CLINVAR_VCF},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    
    &> ${LOG}

if [ $? -eq 0 ]; then
    echo "VEP annotation completed for sample ${SAMPLE}"
else
    echo "VEP annotation failed for sample ${SAMPLE}"
fi
