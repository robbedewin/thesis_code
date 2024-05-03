rule gridss_calling:
    input:
        tumor_bam="results/bam/{sample}_tumor.bam",
        tumor_bai="results/bam/{sample}_tumor.bam.bai",
        normal_bam="results/bam/{sample}_normal.bam",
        normal_bai="results/bam/{sample}_normal.bam.bai",
        reference="resources/genome.fa",
        reference_index="resources/genome.fa.fai",
        blacklist="resources/T2T.excluderanges.bed",
        pon="resources/GRCh38_PoN.vcf",
    output:
        vcf="results/gridss/{sample}_gridss.vcf",
    params:
        working_dir="temp/gridss/{sample}/",
    log:
        "logs/gridss/{sample}.log",
    shell:
        """
        gridss --reference {input.reference} \
               --output {output.vcf} \
               --assembly {params.working_dir}{sample}.gridss.assembly.bam \
               --workingdir {params.working_dir} \
               --blacklist {input.blacklist} \
               --pon {input.pon} \
               --tumor {input.tumor_bam} \
               --normal {input.normal_bam} \
               2> {log}
        """

rule call_variants:
    input:
        bam="results/recal/{sample}/{sample}_dna_normal_recal.bam",
        reference="resources/genome.fa",
    output:
        vcf="results/gridss/vcfs/{sample}.vcf.gz"
    shell:
        """
        # Command for variant calling with GRIDSS
        # This is a placeholder. Replace with your actual GRIDSS variant calling command.
        gridss -r {input.ref_genome} -o {output.vcf} -a {input.bam} --jvmheap 8g
        """

rule generate_pon:
    input:
        vcf_list=expand("results/gridss/vcfs/{sample}.vcf.gz", sample=unique_samples),
        reference="resources/genome.fa"
    output:
        pon_breakpoint="results/gridss/pondir/gridss_pon_breakpoint.bedpe",
        pon_single_breakend="results/gridss/pondir/gridss_pon_single_breakend.bed"
    params:
        pon_dir="pondir"
    conda:
        "envs/gridss-env.yml"
    shell:
        """
        mkdir -p {params.pon_dir}
        java -Xmx8g -cp $(dirname $(which gridss))/../share/gridss-*/gridss.jar gridss.GeneratePonBedpe \
            $(echo {input.vcf_list} | tr ' ' '\n' | awk '{{print "INPUT=" $0}}') \
            O={output.pon_breakpoint} \
            SBO={output.pon_single_breakend} \
            REFERENCE_SEQUENCE={input.reference}
        """
