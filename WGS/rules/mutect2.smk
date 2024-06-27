# Calls variants in normal samples using Mutect2 to create a Panel of Normals (PoN)
rule mutect_pon1:
    input:
        genome="resources/genome.fa",
        normal="results/recal/{sample}/{sample}_dna_normal_recal.bam",
    output:
        vcf="results/mutect2/pon/{sample}_normal_for_pon.vcf.gz",
    threads: 4
    log:
        "logs/mutect2/pon/{sample}.log",    
    shell:
        """
        gatk --java-options "-Xmx200g -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
        Mutect2 \
            -R {input.genome} \
            -I {input.normal} \
            --max-mnp-distance 0 \
            -O {output.vcf} \
            2>{log}
        """

# Create a sample map for the GenomicsDBImport
rule create_sample_map:
    input:
        expand("results/mutect2/pon/{sample}_normal_for_pon.vcf.gz", sample=unique_samples),
    output:
        "results/mutect2/pon/sample_map",
    run:
        unique_samples = units['sample_name'].unique()  # Ensure unique sample names
        with open(output[0], "w") as map_file:
            for sample in unique_samples:
                vcf_path = f"results/mutect2/pon/{sample}_normal_for_pon.vcf.gz"
                map_file.write(f"{sample}\t{vcf_path}\n")

# Uses GenomicsDBImport to integrate the VCFs from multiple normal samples into a single GenomicsDB workspace
rule mutect_pon2:
    input:
        genome="resources/genome.fa",
        sample_map="results/mutect2/pon/sample_map",
        intervals="resources/whole_genome.intervals",
    output:
        directory("results/mutect2/pon/pon_db"),
    log:
        "logs/mutect2/pon/genomicsdbImport.log",
    shell:
        """
        gatk GenomicsDBImport \
            -R {input.genome} \
            -L {input.intervals} \
            --genomicsdb-workspace-path {output} \
            --sample-name-map {input.sample_map} \
            2>{log}
        """

# Executes CreateSomaticPanelOfNormals to consolidate variant calls from the GenomicsDB workspace into a PoN VCF.
# The PoN helps in identifying recurrent technical artifacts, improving the accuracy of somatic variant calling.
rule mutect_pon3:
    input:
        genome="resources/genome.fa",
        pon_db="results/mutect2/pon/pon_db",
    output:
        pon_vcf="results/mutect2/pon/mutect2.pon.vcf.gz",
    params:
        germline="resources/gnomad_chr.vcf.gz",
    log:
        "logs/mutect2/pon/create_somatic_panel_of_normals.log",
    shell:
        """
        gatk CreateSomaticPanelOfNormals \
            -R {input.genome} \
            --germline-resource {params.germline} \
            -V gendb://{input.pon_db} \
            -O {output.pon_vcf} \
            2>{log}
        """

# Primary rule for somatic variant calling in tumor-normal pairs using Mutect2.
# Outputs raw variant calls in VCF format and generates f1r2.tar.gz files for each pair, aiding in further bias analysis.
rule mutect2: 
    input:
        tumor="results/recal/{sample}/{sample}_dna_tumor_recal.bam",
        normal="results/recal/{sample}/{sample}_dna_normal_recal.bam",
        pon="results/mutect2/pon/mutect2.pon.vcf.gz",
    output:
        vcf="results/mutect2/{sample}/{sample}_variants.vcf",
        f1r2="results/mutect2/{sample}/{sample}_f1r2.tar.gz",
    params:
        genome="resources/genome.fa",
        germline_resource="resources/gnomad_chr.vcf.gz",
    log:
        "logs/mutect2/{sample}/{sample}_variants.log",
    threads: 16,
    shell: 
        """
        gatk Mutect2 \
             -R {params.genome} \
             -I {input.tumor} \
             -I  {input.normal} \
             -normal {wildcards.sample}_normal \
             --germline-resource {params.germline_resource} \
             --panel-of-normals {input.pon} \
             -O {output.vcf} \
             --f1r2-tar-gz {output.f1r2} \
             --native-pair-hmm-threads {threads} \
            &> {log}
        """     


# Utilizes the f1r2.tar.gz files from the Mutect2 rule to model read orientation bias using GATK's LearnReadOrientationModel.
# The model generated is used in subsequent analysis to correct for orientation bias, improving variant call accuracy.
rule learn_read_orientation_model:
    input:
        f1r2=expand("results/mutect2/{sample}/{sample}_f1r2.tar.gz", sample=unique_samples),
    output:
        model="results/mutect2/read_orientation_model.tar.gz",
    log:
        "logs/mutect2/learn_read_orientation_model.log",
    shell:
        """
        # Create input files list
        echo {input.f1r2} | tr ' ' '\\n' | awk '{{print "-I "$$0}}' | paste -sd' ' > input_files.txt

        # Execute GATK command
        gatk LearnReadOrientationModel \
             -I $(cat input_files.txt) \
             -O {output.model} \
             &> {log}

        # Clean up temporary file
        rm input_files.txt
        """


# Generates summaries of read pileup information from normal and tumor samples using GATK's GetPileupSummaries.
# This information is vital for estimating cross-sample contamination in the CalculateContamination step.
rule get_pileup_summaries:
    input:
        bam="results/recal/{sample}/{sample}_dna_{alias}_recal.bam",
        intervals="resources/whole_genome.intervals",
        germline_resource="resources/gnomad_chr.vcf.gz",
    output:
        summaries="results/mutect2/{sample}/{sample}_{alias}_pileup_summaries.table",
    log:
        "logs/mutect2/{sample}/{sample}_{alias}_get_pileup_summaries.log",
    shell:
        """
        gatk GetPileupSummaries \
             -I {input.bam} \
             -L {input.intervals} \
             -V {input.germline_resource} \
             -O {output.summaries} \
             &> {log}
        """

# Calculates the level of cross-sample contamination for each sample based on pileup summaries.
# The contamination estimates are used to adjust variant allele fractions, ensuring more accurate variant calling.
rule calculate_contamination:
    input:
        summaries_normal="results/mutect2/{sample}/{sample}_normal_pileup_summaries.table",
        summaries_tumor="results/mutect2/{sample}/{sample}_tumor_pileup_summaries.table",
    output:
        contamination="results/mutect2/{sample}/{sample}_contamination.table",
    log:
        "logs/mutect2/{sample}/{sample}_calculate_contamination.log",
    shell:
        """
        gatk CalculateContamination \
             -I {input.summaries_tumor} \
             -matched {input.summaries_normal} \
             -O {output.contamination} \
             &> {log}
        """

# Applies GATK's FilterMutectCalls to refine the raw variant calls based on various criteria, including statistical thresholds.
# This step enhances the reliability of variant calls by filtering out potential false positives.
rule filter_mutect2:
    input:
        vcf="results/mutect2/{sample}/{sample}_variants.vcf",
        contamination="results/mutect2/{sample}/{sample}_contamination.table",
        model="results/mutect2/read_orientation_model.tar.gz",
    output:
        filtered_vcf="results/mutect2/{sample}/{sample}_filtered_variants.vcf",
    params:
        genome="resources/genome.fa",
    log:
        "logs/mutect2/{sample}/{sample}_filtered.log",
    shell:
        """
        gatk FilterMutectCalls \
            -R {params.genome} \
            -V {input.vcf} \
            --contamination-table {input.contamination} \
            --ob-priors {input.model} \
            -O {output.filtered_vcf} \
            &> {log}
        """


# Selects variants marked as PASS from the filtered VCF files using GATK's SelectVariants.
# This step isolates high-confidence variant calls for downstream analysis.
rule select_pass_variants:
    input:
        vcf="results/mutect2/{sample}/{sample}_filtered_variants.vcf",
    output:
        vcf="results/mutect2/{sample}/{sample}_pass_variants.vcf",
    log:
        "logs/mutect2/{sample}/{sample}_select_pass_variants.log",
    shell:
        """
        gatk SelectVariants \
             -V {input.vcf} \
             --exclude-filtered \
             -O {output.vcf} \
             &> {log}
        """


rule liftover_chm13_to_grch38_vcf:
    input:
        original_vcf="results/mutect2/{sample}/{sample}_pass_variants.vcf"
    output:
        lifted="results/mutect2/{sample}/{sample}_pass_variants_lifted.vcf",
        unlifted="results/mutect2/{sample}/{sample}_pass_variants_unlifted.vcf",
    params:
        genome="resources/GRCh38/GRCh38.fasta", 
        chain="resources/liftover/chm13v2-to-grch38.chain",
    log:
        "logs/mutect2/{sample}/{sample}_liftover.log",
    shell:
        """
        picard -Xmx4g  LiftoverVcf \
            -I {input.original_vcf} \
            -O {output.lifted} \
            -CHAIN {params.chain} \
            -REJECT {output.unlifted} \
            -R {params.genome} \
            --LIFTOVER_MIN_MATCH 0.95 \
            &> {log}
        """

# Annotates the final PASS variants with functional information using GATK's Funcotator.
rule funcotate_variants:
    input:
        vcf="results/mutect2/{sample}/{sample}_pass_variants_lifted.vcf"
    output:
        annotated_vcf="results/mutect2/{sample}/{sample}_annotated_lifted_variants.maf",
    params:
        genome="resources/GRCh38/GRCh38.fasta",
        data_sources_path="/staging/leuven/stg_00096/references/GRCh38.alt-masked-V2/annotation/funcotator/funcotator_dataSources.v1.7.20200521s/",
    log:
        "logs/mutect2/{sample}/{sample}_funcotate_variants.log",
    shell:
        """
        gatk Funcotator \
            -R {params.genome} \
            -V {input.vcf} \
            -O {output.annotated_vcf} \
            --data-sources-path {params.data_sources_path} \
            --ref-version hg38 \
            --output-file-format MAF \
            &> {log}
        """
