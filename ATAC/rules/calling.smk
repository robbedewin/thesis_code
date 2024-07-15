rule bigWig:
    input:
        bam="results/samtools/{sample}/{sample}_atac.bam",
        bai="results/samtools/{sample}/{sample}_atac.bam.csi", #changed to .csi instead of .bai
    output:
        bigWig="results/deepTools/{sample}/{sample}.bigWig",
    log:
        "logs/deepTools/{sample}_bigWig.log",    
    params:
        genome_size = config["genome_size"],
    shell:
        """
        bamCoverage --bam {input.bam} \
            -p max --binSize 10  --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} --extendReads 175 \
            -o {output.bigWig} \
            &> {log}
        """

rule callpeak:
    input:
        bam="results/samtools/{sample}/{sample}_atac.bam",
        bai="results/samtools/{sample}/{sample}_atac.bam.csi",  #changed to .csi instead of .bai
    output:
        peaks="results/macs2/{sample}/{sample}_peaks.narrowPeak",
        summits="results/macs2/{sample}/{sample}_summits.bed",
        xls="results/macs2/{sample}/{sample}_peaks.xls",
    log:
        "logs/macs2/{sample}_callpeak.log",
    params:
        genome_size = config["genome_size"],
    shell:
        """
        macs2 callpeak -t {input.bam} -f BAM -g {params.genome_size} -n {wildcards.sample} \
            --outdir results/macs2/{wildcards.sample} \
            &> {log}
        """