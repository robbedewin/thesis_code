wildcard_constraints:
    read = "R[12]"

rule fastqc_raw:
    input:
        lambda wildcards: get_fastq_dna(wildcards)[wildcards.read]
    output:
        html="results/qc/fastqc/{sample}/{sample}_dna_{alias}_raw_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_dna_{alias}_raw_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_dna_{alias}_raw_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"


# rule trim_reads_for_fastqc:
#     input:
#         unpack(get_fastq_dna),
#     output:
#         trimmed_r1=temp("results/qc/fastp/{sample}/{sample}_dna_{alias}_R1.trimmed.fq"),
#         trimmed_r2=temp("results/qc/fastp/{sample}/{sample}_dna_{alias}_R2.trimmed.fq"),
#         html="logs/qc/fastp/{sample}/{sample}_dna_{alias}.fastp.html",
#         json="logs/qc/fastp/{sample}/{sample}_dna_{alias}.fastp.json",
#     log:
#         fastp_log="logs/fastp/{sample}/{sample}_dna_{alias}.fastp.log"
#     shell:
#         """
#         fastp -w 12 -i {input.R1} -I {input.R2} \
#               -o {output.trimmed_r1} -O {output.trimmed_r2} \
#               -j {output.json} -h {output.html} &> {log.fastp_log}
#         """

rule fastqc_trimmed:
    input:
        "results/fastp/{sample}/{sample}_dna_{alias}_{read}.trimmed.fq"
    output:
        html="results/qc/fastqc/{sample}/{sample}_dna_{alias}_trimmed_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_dna_{alias}_trimmed_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_dna_{alias}_trimmed_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"

rule samtools_idxstats:
    input:
        bam="results/recal/{sample}/{sample}_dna_{alias}_recal.bam",
        idx="results/recal/{sample}/{sample}_dna_{alias}_recal.bai",
    output:
        "results/qc/samtools/{sample}/{sample}_dna_{alias}_recal.bam.idxstats",
    log:
        "logs/samtools/{sample}/{sample}_dna_{alias}_recal_idxstats.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v3.7.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        bam="results/recal/{sample}/{sample}_dna_{alias}_recal.bam",
    output:
        "results/qc/samtools/{sample}/{sample}_dna_{alias}_recal.bam.stats",
    log:
        "logs/samtools/{sample}/{sample}_dna_{alias}_recal.log",
    wrapper:
        "v3.7.0/bio/samtools/stats"

rule multiqc:
    input:
        expand(
            [
                "results/qc/fastqc/{sample}/{sample}_dna_{alias}_raw_{read}_fastqc.zip",
                "results/qc/fastqc/{sample}/{sample}_dna_{alias}_trimmed_{read}_fastqc.zip",
                "results/qc/samtools/{sample}/{sample}_dna_{alias}_recal.bam.idxstats",
                "results/qc/samtools/{sample}/{sample}_dna_{alias}_recal.bam.stats",
                "results/qc/dedup/{sample}/{sample}_dna_{alias}_duplication_metrics.txt",
            ],
            sample=unique_samples,
            alias=unique_aliases,
            read=["R1", "R2"]
        ),
    output:
        report(
            "results/qc/multiqc_34samples.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    params:
        extra="--fn_as_s_name",  #Use log file names as sample names
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.7.0/bio/multiqc"

