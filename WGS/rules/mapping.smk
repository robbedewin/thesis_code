rule trim_reads:
    input:
        unpack(get_fastq_dna),
    output:
        trimmed_r1=temp("results/fastp/{sample}/{sample}_dna_{alias}_R1.trimmed.fq"),
        trimmed_r2=temp("results/fastp/{sample}/{sample}_dna_{alias}_R2.trimmed.fq"),
        html="logs/fastp/{sample}/{sample}_dna_{alias}.fastp.html",
        json="logs/fastp/{sample}/{sample}_dna_{alias}.fastp.json",
    log:
        fastp_log="logs/fastp/{sample}/{sample}_dna_{alias}.fastp.log"
    shell:
        """
        fastp -w 12 -i {input.R1} -I {input.R2} \
              -o {output.trimmed_r1} -O {output.trimmed_r2} \
              -j {output.json} -h {output.html} &> {log.fastp_log}
        """

rule bwa_map:
    input:
        r1="results/fastp/{sample}/{sample}_dna_{alias}_R1.trimmed.fq",
        r2="results/fastp/{sample}/{sample}_dna_{alias}_R2.trimmed.fq",
        genome="resources/genome.fa",
        idx=multiext("resources/genome.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        sam=temp("results/mapped_reads/{sample}/{sample}_dna_{alias}.sam")
    log:
        map_log="logs/mapped_reads/{sample}/{sample}_dna_{alias}.log.bwamem"
    threads: 18
    shell:
        """
        read -r first_line < <(cat {input.r1} | head -n 1)
        IFS=':' read -r -a rgdata <<< "$first_line"
        FLOWCELLID=${{rgdata[2]}}
        LANE=${{rgdata[3]}}
        BC=${{rgdata[-1]}}
        RG="@RG\\tID:${{FLOWCELLID}}.${{LANE}}\\tPL:ILLUMINA\\tPU:${{FLOWCELLID}}.${{LANE}}.${{BC}}\\tLB:LIB-{wildcards.sample}_{wildcards.alias}-1\\tSM:{wildcards.sample}_{wildcards.alias}"

        bwa-mem2 mem -o {output.sam} -t {threads} -R "$RG" \
        {input.genome} {input.r1} {input.r2} &> {log.map_log}
        """

rule samtools_sort:
    input:
        sam="results/mapped_reads/{sample}/{sample}_dna_{alias}.sam",
        genome="resources/genome.fa"
    output:
        bam=temp("results/mapped_reads/{sample}/{sample}_dna_{alias}.sorted.bam"),
    log:
        sort_log="logs/samtools/{sample}/{sample}_dna_{alias}_sort.log"
    params:
        tmp_prefix="sorttmp_{wildcards.sample}_{wildcards.datatype}_{wildcards.alias}"
    shell:
        """
        samtools sort -@ 18 --write-index \
        -m 2G \
        --reference {input.genome} \
        -T {params.tmp_prefix} \
        -o {output.bam} {input.sam} &> {log.sort_log}
        """

rule mark_duplicates:
    input:
        bam="results/mapped_reads/{sample}/{sample}_dna_{alias}.sorted.bam"
    output:
        marked_bam=temp("results/mapped_reads/{sample}/{sample}_dna_{alias}.marked.bam"),
        metrics="results/qc/dedup/{sample}/{sample}_dna_{alias}_duplication_metrics.txt"
    log:
        markdup_log="logs/picard/{sample}/{sample}_dna_{alias}_markdup.log"
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.marked_bam} \
            M={output.metrics} \
            &> {log.markdup_log}
        """

rule base_recal_table:
    input:
        bam="results/mapped_reads/{sample}/{sample}_dna_{alias}.marked.bam",
        genome="resources/genome.fa",
        genome_index="resources/genome.fa.fai",
        genome_dict="resources/genome.dict",
        known_sites="resources/dbSNPv155_common.vcf.gz",
        known_sites_index="resources/dbSNPv155_common.vcf.gz.tbi"
    output:
        recal_table="results/recal/{sample}/{sample}_dna_{alias}_recal_data.table"
    log:
        recal_log="logs/recal/{sample}/{sample}_dna_{alias}_recal.log"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.genome} \
            --known-sites {input.known_sites} \
            -O {output.recal_table} \
            &> {log.recal_log}
        """

rule apply_bqsr:
    input:
        bam="results/mapped_reads/{sample}/{sample}_dna_{alias}.marked.bam",
        recal_table="results/recal/{sample}/{sample}_dna_{alias}_recal_data.table",
        genome="resources/genome.fa"
    output:
        recal_bam=protected("results/recal/{sample}/{sample}_dna_{alias}_recal.bam")
    log:
        apply_recal_log="logs/recal/{sample}/{sample}_dna_{alias}_apply_recal.log"
    shell:
        """
        gatk ApplyBQSR \
            -R {input.genome} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.recal_bam} \
            --create-output-bam-index false \
            &> {log.apply_recal_log}
        """
# added the --create-output-bam-index false flag to prevent the creation of the index file due to this error: htsjdk.samtools.SAMException: Not creating BAM index

rule samtools_index:
    input:
        bam="results/recal/{sample}/{sample}_dna_{alias}_recal.bam"
    output:
        bai="results/recal/{sample}/{sample}_dna_{alias}_recal.bai"
    log:
        index_log="logs/samtools/{sample}/{sample}_dna_{alias}_index.log"
    shell:
        """
        samtools index -@ 18 {input.bam} {output.bai} &> {log.index_log}
        """
