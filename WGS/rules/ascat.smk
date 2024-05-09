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
        qc_output="results/ascat/qc/summary.csv",
        seg_plot="results/ascat/qc/segments_plot.png",
    params:
        sample_ids=unique_samples,
        ascat_dir="results/ascat/"
    log:
        "logs/ascat/qc.log"
    script:
        "scripts/ascat_qc.R"

rule test:
    input:
        "results/recal/P013/P013_dna_tumor_recal.bam"
    output:
        "results/test/output.txt"
    params:
        sample_ids=unique_samples
    script:
        "scripts/test.R"

