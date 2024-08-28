#!/bin/bash

# Define variables
SAMPLE="P020"
VCF_INPUT="results/vep/${SAMPLE}/${SAMPLE}_gridss_high_confidence_somatic_annotated.vcf"
MAF_OUTPUT="results/vcf2maf/${SAMPLE}.maf"
TUMOR_ID="${SAMPLE}_tumor"
NORMAL_ID="${SAMPLE}_normal"
VEP_PATH="/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/VEP/bin"
VEP_DATA="resources/vep/cache"
REF_FASTA="resources/genome.fa"
VCF2MAF_LOG="logs/vcf2maf/${SAMPLE}.log"

# Create directories if they do not exist
mkdir -p $(dirname ${MAF_OUTPUT})
mkdir -p $(dirname ${VCF2MAF_LOG})

# Run vcf2maf to convert annotated VCF to MAF
echo "Converting VCF to MAF..."
perl ${VEP_PATH}/vcf2maf.pl \
    --input-vcf ${VCF_INPUT} \
    --output-maf ${MAF_OUTPUT} \
    --tumor-id ${TUMOR_ID} \
    --normal-id ${NORMAL_ID} \
    --vcf-tumor-id ${TUMOR_ID} \
    --vcf-normal-id ${NORMAL_ID} \
    --ref-fasta ${REF_FASTA} \
    --ncbi-build T2T-CHM13v2.0 \
    --retain-info gnomAD_AF,gnomAD_AC,gnomAD_AN,gnomAD_AF_afr,gnomAD_AF_ami,gnomAD_AF_amr,gnomAD_AF_asj,gnomAD_AF_eas,gnomAD_AF_fin,gnomAD_AF_nfe,gnomAD_AF_oth,gnomAD_AF_sas,dbSNP_RS,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN \
    --inhibit-vep \
    > ${VCF2MAF_LOG} 2>&1

echo "VCF to MAF conversion completed for sample ${SAMPLE}"
