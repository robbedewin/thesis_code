rule aggregate_counts:
    input:
        expand("results/star/{sample}/{sample}_ReadsPerGene.out.tab", sample=units["sample_name"].unique())
    output:
        counts="results/counts/combined_counts.tsv"
    log:
        "logs/aggregate_counts.log"
    params:
        sample_info=sample_info,
    script:
        "scripts/aggregate_counts.py"

rule deseq2:
    input:
        counts="results/counts/combined_counts.tsv",
    output:
        results_dds="results/deseq2/deseq2_results.tsv",
    params:
        sample_info_file="results/sample_info.csv",
        plot_dir="results/deseq2/plots/",
    log:
        "logs/deseq2.log"
    script:
        "scripts/deseq2.R"


    