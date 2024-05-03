rule trim_atac_reads:
    input:
        unpack(get_fastq_atac) 
    output:
        trimmed_r1="results/fastp/{sample}/{sample}_atac_R1.trimmed.fq",
        trimmed_r2="results/fastp/{sample}/{sample}_atac_R2.trimmed.fq",
        html="logs/fastp/{sample}/{sample}_atac.fastp.html",
        json="logs/fastp/{sample}/{sample}_atac.fastp.json",
    log:
        fastp_log="logs/fastp/{sample}/{sample}_atac.fastp.log"
    shell:
        """
        fastp -w 12 -i {input.R1} -I {input.R2} \
              -o {output.trimmed_r1} -O {output.trimmed_r2} \
              -j {output.json} -h {output.html} &> {log.fastp_log}
        """

rule minimap2:
    input:
        r1="results/fastp/{sample}/{sample}_atac_R1.trimmed.fq",
        r2="results/fastp/{sample}/{sample}_atac_R2.trimmed.fq",
        reference_index="resources/genome.mmi"
    output:
        sam="results/minimap2/{sample}/{sample}_aligned.sam"
    params:
        preset="-ax sr", #'sr' for short reads
        extra="", 
    threads: 8 
    log:
        "logs/minimap2/{sample}_minimap2.log"
    shell:
        """
        minimap2 {params.preset} {params.extra} \
                 -t {threads} \
                 {input.reference_index} {input.r1} {input.r2} > {output.sam} \
                 2> {log}
        """

rule bwa_mem2:
    input:
        r1="results/fastp/{sample}/{sample}_atac_R1.trimmed.fq",
        r2="results/fastp/{sample}/{sample}_atac_R2.trimmed.fq",
        ref_genome="resources/genome.fa",
        idx=multiext("resources/genome.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        sam="results/bwamem2/{sample}/{sample}_aligned.sam"
    params:
        extra="",  # Placeholder for any additional BWA-MEM2 parameters you might need
    threads: 18  # Adjust based on your system's resources
    log:
        "logs/bwamem2/{sample}_bwamem2.log"
    shell:
        """
        read -r first_line < <(cat {input.r1} | head -n 1)
        IFS=':' read -r -a rgdata <<< "$first_line"
        FLOWCELLID=${{rgdata[2]}}
        LANE=${{rgdata[3]}}
        BC=${{rgdata[-1]}}
        RG="@RG\\tID:${{FLOWCELLID}}.${{LANE}}\\tPL:ILLUMINA\\tPU:${{FLOWCELLID}}.${{LANE}}.${{BC}}\\tLB:LIB-{wildcards.sample}_atac-1\\tSM:{wildcards.sample}_atac"

        bwa-mem2 mem -o {output.sam} -t {threads} -R "$RG" \
        {input.ref_genome} {input.r1} {input.r2} &> {log}
        """


rule samtools_sort_index:
    input:
        sam="results/bwamem2/{sample}/{sample}_aligned.sam",
        genome="resources/genome.fa",
    output:
        bam="results/samtools/{sample}/{sample}_atac.bam",
        csi="results/samtools/{sample}/{sample}_atac.bam.csi",
    log:
        sort_log="logs/samtools/{sample}/{sample}_sort.log",
        index_log="logs/samtools/{sample}/{sample}_index.log",
    params:
        tmp_prefix="sorttmp_{wildcards.sample}",
    threads: 18
    shell:
        """
        samtools sort --threads {threads} --write-index -m 2G \
        --reference {input.genome} -T {params.tmp_prefix} \
        -o {output.bam} {input.sam} &> {log.sort_log}

        samtools index --threads {threads} \
        {output.bam} --csi {output.csi} &> {log.index_log}
        """