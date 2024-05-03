wildcard_constraints:
    read = "R[12]"

rule fastqc_raw:
    input:
        lambda wildcards: get_fastq_rna(wildcards)[wildcards.read]
    output:
        html="results/qc/fastqc/{sample}/{sample}_rna_raw_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_rna_raw_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_rna_raw_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"


rule fastqc_trimmed:
    input:
        "results/fastp/{sample}/{sample}_rna_{read}.trimmed.fq"
    output:
        html="results/qc/fastqc/{sample}/{sample}_rna_trimmed_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_rna_trimmed_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_rna_trimmed_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"

rule multiqc:
    input:
        expand(
            [
                "results/qc/fastqc/{sample}/{sample}_rna_raw_{read}_fastqc.zip",
                "results/qc/fastqc/{sample}/{sample}_rna_trimmed_{read}_fastqc.zip",
            ],
            sample=unique_samples,
            read=["R1", "R2"]
        ),
    output:
        report(
            "results/qc/multiqc_RNA.html",
            caption="../report/multiqc.rst",
            category="Quality control RNA",
        ),
    params:
        extra="--fn_as_s_name",  #Use log file names as sample names
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.7.0/bio/multiqc"

