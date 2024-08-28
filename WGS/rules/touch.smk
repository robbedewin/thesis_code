# Ensure touch_old_files is run before any other rule
ruleorder: touch_old_files > all

old_samples = ["P013", "P016", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P028"]
aliases = ["tumor", "normal"]

# Define the outputs that need to be marked as completed
files_to_touch = [
    # Mapping with bwa mem
    expand("results/mapped_reads/{sample}/{sample}_dna_{alias}.sam", sample=old_samples, alias=aliases),
    # Deduplication
    expand("results/mapped_reads/{sample}/{sample}_dna_{alias}.marked.bam", sample=old_samples, alias=aliases),
    # Recalibration (table and applying)
    expand("results/recal/{sample}/{sample}_dna_{alias}_recal.bam", sample=old_samples, alias=aliases),
    # MultiQC report (fastqc, samtools_(idx)stats, duplication_metrics)
    "results/qc/multiqc_test.html",
    # Ascat
    expand("results/ascat/{sample}/", sample=old_samples),
    # Mutect2 PON
    "results/mutect2/pon/mutect2.pon.vcf.gz",
    # Mutect2
    expand("results/mutect2/{sample}/{sample}_variants.vcf", sample=old_samples),
    # Mutect2 filtering and PASS variants
    expand("results/mutect2/{sample}/{sample}_pass_variants.vcf", sample=old_samples),
    # Mutect2 liftover and annotation
    expand("results/mutect2/{sample}/{sample}_pass_variants_lifted.vcf", sample=old_samples),
    expand("results/mutect2/{sample}/{sample}_annotated_lifted_variants.maf", sample=old_samples),
    # Svaba
    expand("results/svaba/{sample}/{sample}_svaba.somatic.sv.vcf", sample=old_samples),
    # Gridss PON
    "results/gridss/pondir/gridss_pon_breakpoint.bedpe",
    # Gridss SV calling
    expand("results/gridss/{sample}/{sample}_all_calls.vcf", sample=old_samples),
    # Gridss filtering
    expand("results/gridss/{sample}/{sample}_high_confidence_somatic.vcf.bgz", sample=old_samples)
]

recal_files = [
    expand("results/recal/{sample}/{sample}_dna_{alias}_recal.bam", sample=old_samples, alias=aliases),
    
]


rule touch_old_files:
    output:
        touch_files=recal_files
    run:
        for file in output.touch_files:
            shell(f"touch {file}")




# Marking files as ancient for old samples
# for sample in old_samples:
#     rule:
#         output:
#             ancient(f"results/mapped_reads/{sample}/{sample}_dna_tumor.sam"),
#             ancient(f"results/mapped_reads/{sample}/{sample}_dna_normal.sam"),
#             ancient(f"results/mapped_reads/{sample}/{sample}_dna_tumor.marked.bam"),
#             ancient(f"results/mapped_reads/{sample}/{sample}_dna_normal.marked.bam")