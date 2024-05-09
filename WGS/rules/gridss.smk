# GRIDSS processing for normal samples
rule gridss_normal:
    input:
        bam = "results/recal/{sample}/{sample}_dna_normal_recal.bam",
        blacklist = "resources/T2T.excluderanges.bed"
    output:
        vcf = "results/gridss/{sample}_normal.vcf",
        assembly = "results/gridss/{sample}_normal.assembly.bam"
    params:
        ref = "resources/genome.fa",
        jar = "/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar",
        working_dir = "results/gridss/temp/{sample}/",
        tmp_dir = "results/gridss/temp/{sample}/tmp"
    threads: 8
    log: 
        "logs/gridss/{sample}_normal.log"
    shell:
        """
        mkdir -p {params.working_dir} {params.tmp_dir}
        java -Xmx30g -Dsamjdk.create_index=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true \
        -Djava.io.tmpdir={params.tmp_dir} \
        -jar {params.jar} \
        --reference {params.ref} \
        --output {output.vcf} \
        --assembly {output.assembly} \
        --blacklist {input.blacklist} \
        --threads {threads} \
        --workingdir {params.working_dir} \
        {input.bam}
        """

# Generation of Panel of Normals (PON)
rule generate_pon:
    input:
        vcfs=expand("results/gridss/{sample}_normal.vcf", sample=unique_samples),
        
    output:
        pon_breakpoint="results/gridss/pondir/gridss_pon_breakpoint.bedpe",
        pon_single_breakend="results/gridss/pondir/gridss_pon_single_breakend.bed"
    params:
        jar="/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar",
        pon_dir="results/gridss/pondir",
        ref="resources/genome.fa",
    shell:
        """
        mkdir -p {params.pon_dir}
        java -Xmx8g \
            -cp {params.jar} \
            gridss.GeneratePonBedpe \
            $(echo {input.vcfs} | tr ' ' '\n' | awk '{{print "INPUT=" $0}}') \
            O={output.pon_breakpoint} \
            SBO={output.pon_single_breakend} \
            REFERENCE_SEQUENCE={params.ref}
        """

# Somatic structural variant calling
rule call_somatic_structural_variants:
    input:
        normal="results/recal/{sample}/{sample}_dna_normal_recal.bam",
        tumor="results/recal/{sample}/{sample}_dna_tumor_recal.bam",
    output:
        vcf="results/gridss/{sample}/{sample}_all_calls.vcf"
    params:
        ref="resources/genome.fa",
        blacklist="resources/T2T.excluderanges.bed",
        jar="/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar"
    threads: 8
    log:
        small_log="logs/gridss/{sample}_somatic_SV.log",
        full_log="logs/gridss/{sample}_somatic_SV.full.log"
    shell:
        """     
        gridss \
            -r {params.ref} \
            -j  {params.jar} \
            -o {output.vcf} \
            -b {params.blacklist} \
            --skipsoftcliprealignment \
            {input.normal} \
            {input.tumor} \
            &> {log.small_log}
            mv gridss.full.*.log {log.full}
        """

rule gridss_somatic_filter:
    input:
        all_calls="results/gridss/{sample}/{sample}_all_calls.vcf",
        pon_dir="results/gridss/pondir/",
    output:
        high_confidence_somatic="results/gridss/{sample}/{sample}_high_confidence_somatic.vcf.gz",
        high_and_low_confidence_somatic="results/gridss/{sample}/{sample}_high_and_low_confidence_somatic.vcf.gz"
    params:
        scriptdir="/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/"
    log:
        "logs/gridss/{sample}_somatic_filter.log"
    shell:
        """
        gridss_somatic_filter \
            --pondir {input.pon_dir} \
            --input {input.all_calls} \
            --output {output.high_confidence_somatic} \
            --fulloutput {output.high_and_low_confidence_somatic} \
            --scriptdir {params.scriptdir} \
            -n 1 \
            -t 2 \
            &> {log}
        """
