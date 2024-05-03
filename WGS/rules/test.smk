The three steps to create a panel of normals are:
Step 1: Run Mutect2 in tumor-only mode for each normal sample:

    gatk Mutect2 -R reference.fasta -I normal1.bam --max-mnp-distance 0 -O normal1.vcf.gz

    gatk Mutect2 -R reference.fasta -I normal2.bam --max-mnp-distance 0 -O normal2.vcf.gz
... etc.

Step 2: Create a GenomicsDB from the normal Mutect2 calls:

    gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list \
        --genomicsdb-workspace-path pon_db \
        -V normal1.vcf.gz \
        -V normal2.vcf.gz \
        -V normal3.vcf.gz
Step 3: Combine the normal calls using CreateSomaticPanelOfNormals:

    gatk CreateSomaticPanelOfNormals -R reference.fasta \
      --germline-resource af-only-gnomad.vcf.gz \
      -V gendb://pon_db \
      -O pon.vcf.gz

# Define your samples somewhere in the Snakefile
SAMPLES = ["sample1", "sample2", "sample3"]  # Add your sample names here

# Run Mutect2 in tumor-only mode for each normal sample
rule panel_of_normals1:
    input:
        genome="resources/genome.fa",
        normal=lambda wildcards: f"results/recal/{wildcards.sample}/{wildcards.sample}_dna_normal_recal.bam",
    output:
        vcf="results/mutect2/pon/{sample}_normal_for_pon.vcf.gz",
    params:
        max_mnp_distance=0,
    log:
        "logs/mutect2/pon/{sample}.log",
    shell:
        """
        gatk Mutect2 \
            -R {input.genome} \
            -I {input.normal} \
            --max-mnp-distance {params.max_mnp_distance} \
            -O {output.vcf} \
            2>{log}
        """

# Create a GenomicsDB from the normal Mutect2 calls
rule panel_of_normals2:
    input:
        genome="resources/genome.fa",
        vcfs=expand("results/mutect2/pon/{sample}_normal_for_pon.vcf.gz", sample=SAMPLES),
        intervals="resources/intervals.interval_list",
    output:
        pon_db="resources/pon_db",  # This is a directory path for the GenomicsDB workspace
    log:
        "logs/genomicsdbimport.log",
    shell:
        """
        gatk GenomicsDBImport \
            -R {input.genome} \
            -L {input.intervals} \
            --genomicsdb-workspace-path {output.pon_db} \
            {' '.join('-V ' + path for path in {input.vcfs})} \
            2>{log}
        """

# Combine the normal calls using CreateSomaticPanelOfNormals
rule panel_of_normals3:
    input:
        genome="resources/genome.fa",
        germline="resources/af-only-gnomad.vcf.gz",
        pon_db="resources/pon_db",  # Reference the GenomicsDB workspace created in the previous step
    output:
        pon_vcf="results/mutect2/pon/mutect2.pon.vcf.gz",
    log:
        "logs/create_somatic_panel_of_normals.log",
    shell:
        """
        gatk CreateSomaticPanelOfNormals \
            -R {input.genome} \
            --germline-resource {input.germline} \
            -V gendb://{input.pon_db} \
            -O {output.pon_vcf} \
            2>{log}
        """
