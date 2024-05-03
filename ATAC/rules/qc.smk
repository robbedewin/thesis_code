wildcard_constraints:
    read = "R[12]"

rule fastqc_raw:
    input:
        lambda wildcards: get_fastq_atac(wildcards)[wildcards.read]
    output:
        html="results/qc/fastqc/{sample}/{sample}_atac_raw_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_atac_raw_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_atac_raw_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"


rule fastqc_trimmed:
    input:
        "results/fastp/{sample}/{sample}_atac_{read}.trimmed.fq"
    output:
        html="results/qc/fastqc/{sample}/{sample}_atac_trimmed_{read}.html",
        zip="results/qc/fastqc/{sample}/{sample}_atac_trimmed_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}/{sample}_atac_trimmed_{read}.log"
    threads: 1
    resources:
        mem_mb = 10000
    wrapper:
        "v3.7.0/bio/fastqc"


rule samtools_idxstats:
    input:
        bam="results/samtools/{sample}/{sample}_atac.bam",
        idx="results/samtools/{sample}/{sample}_atac.bam.bai",
    output:
        "results/qc/samtools/{sample}/{sample}_atac.bam.idxstats",
    log:
        "logs/samtools/{sample}/{sample}_idxstats.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v3.7.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        bam="results/samtools/{sample}/{sample}_atac.bam",
    output:
        "results/qc/samtools/{sample}/{sample}_atac.bam.stats",
    log:
        "logs/samtools/{sample}/{sample}_stats.log",
    wrapper:
        "v3.7.0/bio/samtools/stats"


rule samtools_flagstat:
    input:
        bam="results/samtools/{sample}/{sample}_atac.bam",
    output:
        flagstat="results/qc/samtools/{sample}/{sample}_atac.flagstat.txt"
    log:
        "logs/samtools/{sample}/{sample}_flagstat.log",
    wrapper:
        "v3.7.0/bio/samtools/flagstat"

rule multiqc:
    input:
        expand(
            [
                "results/qc/fastqc/{sample}/{sample}_atac_raw_{read}_fastqc.zip",
                "results/qc/fastqc/{sample}/{sample}_atac_trimmed_{read}_fastqc.zip",
                "results/qc/samtools/{sample}/{sample}_atac.bam.idxstats",
                "results/qc/samtools/{sample}/{sample}_atac.bam.stats",
                "results/qc/samtools/{sample}/{sample}_atac.flagstat.txt",
            ],
            sample=unique_samples,
            read=["R1", "R2"]
        ),
    output:
        report(
            "results/qc/multiqc_ATAC.html",
            caption="../report/multiqc.rst",
            category="Quality control ATAC",
        ),
    params:
        extra="--fn_as_s_name",  #Use log file names as sample names
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.7.0/bio/multiqc"


