rule ascat:
    input:
        tumor_bam="results/recal/{sample}/{sample}_dna_tumor_recal.bam",
        normal_bam="results/recal/{sample}/{sample}_dna_normal_recal.bam",
    output:
        directory("results/ascat/{sample}/")
    threads: 8
    log:
        "logs/ascat/{sample}.log"
    script:
        "scripts/ascat.R"

rule ascat_qc:
    input:
        expand(["results/ascat/{sample}/"], sample=[u.sample_name for u in units.itertuples()])
    output:
        qc_output="results/ascat/{sample}/{sample}_qc_summary.csv",
        seg_plot="results/ascat/{sample}/{sample}_segments_plot.png"
    log:
        "logs/ascat_processing/{sample}.log"
    script:
        "scripts/process_ascat_output.R"

