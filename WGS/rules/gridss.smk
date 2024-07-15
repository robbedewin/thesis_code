# GRIDSS processing for normal samples
rule gridss_normal:
    input:
        bam = "results/recal/{sample}/{sample}_dna_normal_recal.bam",
        blacklist = "resources/T2T.excluderanges.bed",
    output:
        vcf = "results/gridss/normal/{sample}_normal.vcf",
        assembly = "results/gridss/normal/{sample}_normal.assembly.bam"
    params:
        ref = "resources/genome.fa",
        jar = "/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar",
        working_dir = "results/gridss/temp/{sample}/",
        tmp_dir = "results/gridss/temp/{sample}/tmp",
    threads: 8
    log: 
        "/staging/leuven/stg_00096/home/rdewin/WGS/logs/gridss/normal/{sample}_normal.log"
    shell:
        """
        gridss \
            --reference {params.ref} \
            --jar {params.jar} \
            --output {output.vcf} \
            --assembly {output.assembly} \
            --blacklist {input.blacklist} \
            --threads {threads} \
            --workingdir {params.working_dir} \
            {input.bam} \
            &> {log}
        """

# Generation of Panel of Normals (PON)
rule generate_pon:
    input:
        vcfs=expand("results/gridss/normal/{sample}_normal.vcf", sample=unique_samples),
    output:
        pon_breakpoint="results/gridss/pondir/gridss_pon_breakpoint.bedpe",
        pon_single_breakend="results/gridss/pondir/gridss_pon_single_breakend.bed"
    params:
        jar="/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar",
        pon_dir="results/gridss/pondir",
        ref="resources/genome.fa",
    log:
        "logs/gridss/PON.log"
    shell:
        """
        mkdir -p {params.pon_dir}
        java -Xmx8g \
            -cp {params.jar} \
            gridss.GeneratePonBedpe \
            $(echo {input.vcfs} | tr ' ' '\n' | awk '{{print "INPUT=" $0}}') \
            O={output.pon_breakpoint} \
            SBO={output.pon_single_breakend} \
            NORMAL_ORDINAL=0 \
            REFERENCE_SEQUENCE={params.ref} \
            &> {log}
        """

# Somatic structural variant calling
rule gridss_call_somatic_SV:
    input:
        normal="results/recal/{sample}/{sample}_dna_normal_recal.bam",
        tumor="results/recal/{sample}/{sample}_dna_tumor_recal.bam",
    output:
        vcf="results/gridss/{sample}/{sample}_all_calls.vcf",
    params:
        ref="resources/genome.fa",
        blacklist="resources/T2T.excluderanges.bed",
        working_dir="results/gridss/wrk/{sample}/",
        jar="/lustre1/project/stg_00096/home/rdewin/system/miniconda/envs/WGS/share/gridss-2.13.2-3/gridss.jar"
    threads: 8
    log:
        "logs/gridss/{sample}/{sample}_somatic_SV.log",
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
            --workingdir {params.working_dir} \
            &> {log}
        """

rule gridss_somatic_filter:
    input:
        all_calls="results/gridss/{sample}/{sample}_all_calls.vcf",
        pon_breakpoint="results/gridss/pondir/gridss_pon_breakpoint.bedpe",
    output:
        high_confidence_somatic="results/gridss/{sample}/{sample}_high_confidence_somatic.vcf.bgz",
        high_and_low_confidence_somatic="results/gridss/{sample}/{sample}_high_and_low_confidence_somatic.vcf.bgz",
    params:
        plot_dir="results/gridss/{sample}/plots/",
        pon_dir="results/gridss/pondir/",
        ref="BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
    log:
        "logs/gridss/{sample}/{sample}_somatic_filter.log"
    shell:
        """
        gridss_somatic_filter \
            --pondir {params.pon_dir} \
            --ref {params.ref} \
            --input {input.all_calls} \
            --output {output.high_confidence_somatic} \
            --fulloutput {output.high_and_low_confidence_somatic} \
            --plotdir {params.plot_dir} \
            -n 1 \
            -t 2 \
            &> {log}
        mv {output.high_confidence_somatic}.bgz {output.high_confidence_somatic}
        mv {output.high_confidence_somatic}.bgz.tbi {output.high_confidence_somatic}.tbi
        mv {output.high_and_low_confidence_somatic}.bgz {output.high_and_low_confidence_somatic}
        mv {output.high_and_low_confidence_somatic}.bgz.tbi {output.high_and_low_confidence_somatic}.tbi
        """

rule unzip_bgz_files:
    input:
        high_confidence_somatic="results/gridss/{sample}/{sample}_high_confidence_somatic.vcf.bgz",
        high_and_low_confidence_somatic="results/gridss/{sample}/{sample}_high_and_low_confidence_somatic.vcf.bgz",
    output:
        high_confidence_somatic="results/gridss/{sample}/{sample}_high_confidence_somatic.vcf",
        high_and_low_confidence_somatic="results/gridss/{sample}/{sample}_high_and_low_confidence_somatic.vcf",
    shell:
        """
        bgzip -d -c {input.high_confidence_somatic} > {output.high_confidence_somatic}
        bgzip -d -c {input.high_and_low_confidence_somatic} > {output.high_and_low_confidence_somatic}
        """
