rule trim_reads:
    input:
        unpack(get_fastq)
    output:
        trimmed_r1=temp("results/fastp/{sample}/{unit}_R1.trimmed.fq"),
        trimmed_r2=temp("results/fastp/{sample}/{unit}_R2.trimmed.fq"),
        html="logs/fastp/{sample}/{unit}.fastp.html"
    log:
        fastp_log="logs/fastp/{sample}/{unit}.fastp.log"
    shell:
        """
        fastp -w 12 -i {input.r1} -I {input.r2} \
              -o {output.trimmed_r1} -O {output.trimmed_r2} \
              -h {output.html} &> {log.fastp_log}
        """
    
rule bwa_map:
    input:
        r1="results/fastp/{sample}/{unit}_R1.trimmed.fq",
        r2="results/fastp/{sample}/{unit}_R2.trimmed.fq"
    output:
        sam=temp("results/mapped_reads/{sample}/{unit}.sam")
    log:
        map_log="logs/mapped_reads/{sample}/{unit}.log.bwamem"
    params:
        ref_genome=config['ref_genome_bwa'],
    shell:
        """
        read -r first_line < <(cat {input.r1} | head -n 1)
        IFS=':' read -r -a rgdata <<< "$first_line"
        FLOWCELLID=${{rgdata[2]}}
        LANE=${{rgdata[3]}}
        BC=${{rgdata[-1]}}
        RG="@RG\\tID:${{FLOWCELLID}}.${{LANE}}\\tPL:ILLUMINA\\tPU:${{FLOWCELLID}}.${{LANE}}.${{BC}}\\tLB:LIB-{wildcards.sample}_{wildcards.unit}-1\\tSM:{wildcards.sample}"

        bwa-mem2 mem -o {output.sam} -t 18 -R "$RG" \
        {params.ref_genome} {input.r1} {input.r2} &> {log.map_log}
        """

rule samtools_sort:
    input:
        sam="results/mapped_reads/{sample}/{unit}.sam"
    output:
        bam=temp("results/mapped_reads/{sample}/{unit}.sorted.bam"),
    log:
        sort_log="logs/samtools/{sample}/{unit}_sort.log"
    params:
        ref_genome=config['ref_genome'],
        tmp_prefix="sorttmp_{wildcards.sample}_{wildcards.unit}_sorted"
    shell:
        """
        samtools sort -@ 18 --write-index \
        -m 2G \
        --reference {params.ref_genome} \
        -T {params.tmp_prefix} \
        -o {output.bam} {input.sam} &> {log.sort_log}
        """

rule base_recal_table:
    input:
        bam="results/mapped_reads/{sample}/{unit}.sorted.bam",
    output:
        recal_table="results/recal/{sample}/{unit}_recal_data.table"
    params:
        ref_genome=config['ref_genome'],
        known_sites=config['known_sites']
    log:
        recal_log="logs/recal/{sample}/{unit}_recal.log"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {params.ref_genome} \
            --known-sites {params.known_sites} \
            -O {output.recal_table} \
            &> {log.recal_log}
        """

rule apply_bqsr:
    input:
        bam="results/mapped_reads/{sample}/{unit}.sorted.bam",
        recal_table="results/recal/{sample}/{unit}_recal_data.table"
    output:
        recal_bam="results/recal/{sample}/{unit}_recal.bam"
    params:
        ref_genome=config['ref_genome']
    log:
        apply_recal_log="logs/recal/{sample}/{unit}_apply_recal.log"
    shell:
        """
        gatk ApplyBQSR \
            -R {params.ref_genome} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.recal_bam} \
            &> {log.apply_recal_log}
        """

rule samtools_index:
    input:
        bam="results/recal/{sample}/{unit}_recal.bam"
    output:
        bai=protected("results/recal/{sample}/{unit}_recal.bai")
    log:
        index_log="logs/samtools/{sample}/{unit}_index.log"
    shell:
        """
        samtools index -@ 18 {input.bam} {output.bai} &> {log.index_log}
        """

# rule haplotype_caller:
#     input:
#         bam="results/recal/{sample}/{unit}_recal.bam"
#     output:
#         vcf="results/variants/{sample}/{unit}.vcf"
#     params:
#         ref_genome=config['ref_genome'],
#         known_sites=config['known_sites']
#     log:
#         call_log="logs/variants/{sample}/{unit}_call.log"
#     shell:
#         """
#         gatk HaplotypeCaller \
#             -R {params.ref_genome} \
#             -I {input.bam} \
#             -O {output.vcf} \
#             --dbsnp {params.known_sites} \
#             &> {log.call_log}
#         """

# rule ascat:
#     input:
        
#     output:
#         ascat="results/ascat/{sample}/{unit}.ascat"
#     params:
#         ref_genome=config['ref_genome']
#     log:
#         ascat_log="logs/ascat/{sample}/{unit}_ascat.log"
#     shell:
#         """
#         ascat.R \
#             -v {input.vcf} \
#             -r {params.ref_genome} \
#             -o {output.ascat} \
#             &> {log.ascat_log}
#         """