rule trim_rna_reads:
    input:
        unpack(get_fastq_rna) 
    output:
        trimmed_r1="results/fastp/{sample}/{sample}_rna_R1.trimmed.fq",
        trimmed_r2="results/fastp/{sample}/{sample}_rna_R2.trimmed.fq",
        html="logs/fastp/{sample}/{sample}_rna.fastp.html",
        json="logs/fastp/{sample}/{sample}_rna.fastp.json",
    log:
        fastp_log="logs/fastp/{sample}/{sample}_rna.fastp.log"
    shell:
        """
        fastp -w 12 -i {input.R1} -I {input.R2} \
              -o {output.trimmed_r1} -O {output.trimmed_r2} \
              -j {output.json} -h {output.html} &> {log.fastp_log}
        """

rule STAR_align:
    input:
        r1="results/fastp/{sample}/{sample}_rna_R1.trimmed.fq",
        r2="results/fastp/{sample}/{sample}_rna_R2.trimmed.fq",
        genome_dir="resources/star_genome",
    output:
        bam="results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        quant="results/star/{sample}/{sample}_ReadsPerGene.out.tab",
    threads: 14
    params:
        out_prefix="results/star/{sample}/{sample}_",
    log:
        "logs/star/{sample}/{sample}_align.log"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.genome_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes NH HI AS NM MD \
             --outSAMstrandField intronMotif \
             --quantMode GeneCounts \
             --twopassMode Basic \
        > {log} 2>&1
        """

rule STAR_Fusion:
    input:
        r1="results/fastp/{sample}/{sample}_rna_R1.trimmed.fq",
        r2="results/fastp/{sample}/{sample}_rna_R2.trimmed.fq"
    output:
        dir=directory("results/star_fusion/{sample}"),
        abridged="results/star_fusion/{sample}/star-fusion.fusion_predictions.abridged.tsv",
    threads: 14
    params:
        genome_lib_dir="/staging/leuven/stg_00096/references/CTAT_Resources/T2T-CHM13_CTAT_lib_Feb162023.plug-n-play/ctat_genome_lib_build_dir",
        singularity_image="/staging/leuven/stg_00096/home/rdewin/singularity/star-fusion.v1.12.0.simg",
    log:
        "logs/star_fusion/{sample}/{sample}_fusion.log"
    shell:
        """
        mkdir -p {output.dir}
        singularity exec -e -B `pwd` -B {params.genome_lib_dir} -B /staging/leuven/stg_00096/home/rdewin/RNA/{output.dir} \
        {params.singularity_image} \
         STAR-Fusion \
            --left_fq {input.r1} \
            --right_fq {input.r2} \
            --genome_lib_dir {params.genome_lib_dir} \
            --CPU {threads} \
            --output_dir {output.dir} \
            --FusionInspector validate \
            --examine_coding_effect \
            --denovo_reconstruct \
        > {log} 2>&1
        """


rule samtools_index_RNA:
    input:
        "results/{sample}/{sample}.bam"
    output:
        "results/{sample}/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}/{sample}.log"
    threads: 14
    shell:
        """
        samtools index -@ {threads} {input} {output} &> {log}
        """
