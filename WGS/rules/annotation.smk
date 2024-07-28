# rule get_vep_cache:
#     output:
#         directory("resources/vep/cache"),
#     params:
#         species="homo_sapiens",
#         build="T2T-CHM13v2.0",
#         release="112",
#         url="https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/indexed_vep_cache/Homo_sapiens-GCA_009914755.4-2022_10.tar.gz",
#     log:
#         "logs/vep/cache.log",
#     cache: "omit-software",  # save space and time with between workflow caching (see docs)
#     wrapper:
#         "v3.13.6/bio/vep/cache"


rule download_vep_cache:
    output:
        cache_dir=directory("resources/vep/cache")
    params:
        url="https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_08/indexed_vep_cache/homo_sapiens_gca009914755v4_T2T-CHM13v2.0.tar.gz",
        assembly="T2T-CHM13v2.0"
    log:
        "logs/vep/cache_download.log"
    shell:
        """
        mkdir -p {output.cache_dir}
        wget {params.url} -P {output.cache_dir}
        tar -xzf {output.cache_dir}/homo_sapiens_gca009914755v4_T2T-CHM13v2.0.tar.gz -C {output.cache_dir}
        """

rule annotate_variants_mutect:
    input:
        calls="results/mutect2/{sample}/{sample}_pass_variants_new.vcf",
        cache="resources/vep/cache",
        fasta="resources/genome.fa",
        fai="resources/genome.fa.fai",
        gnomad_vcf="resources/gnomad.vcf.gz",
        dbsnp_vcf="resources/dbSNPv155_common.vcf.gz",
        clinvar_vcf="resources/clinvar.vcf.gz"
    output:
        calls="results/vep/{sample}/{sample}_mutect_pass_variants_annotated.vcf",
        stats="results/vep/{sample}/{sample}_mutect_pass_variants_annotated.html"
    params:
        species="homo_sapiens_gca009914755v4",
        assembly="T2T-CHM13v2.0",
        cache_version=106
    log:
        "logs/vep/{sample}/mutect_pass_variants_annotated.log"
    threads: 8
    shell:
        """
        ./rules/scripts/run_vep_annotation.sh {wildcards.sample} {input.calls} {output.calls} {output.stats} {input.cache} {input.fasta} {input.gnomad_vcf} {input.dbsnp_vcf} {input.clinvar_vcf} {log} {params.cache_version} {params.species} {params.assembly}
        """

rule annotate_variants_gridss:
    input:
        calls="results/gridss/{sample}/{sample}_high_confidence_somatic.vcf",
        cache="resources/vep/cache",
        fasta="resources/genome.fa",
        fai="resources/genome.fa.fai",
        gnomad_vcf="resources/gnomad.vcf.gz",
        dbsnp_vcf="resources/dbSNPv155_common.vcf.gz",
        clinvar_vcf="resources/clinvar.vcf.gz"
    output:
        calls="results/vep/{sample}/{sample}_gridss_high_confidence_somatic_annotated.vcf",
        stats="results/vep/{sample}/{sample}_gridss_high_confidence_somatic_annotated.html"
    params:
        species="homo_sapiens_gca009914755v4",
        assembly="T2T-CHM13v2.0",
        cache_version=106
    log:
        "logs/vep/{sample}/gridss_high_confidence_somatic_annotated.log"
    threads: 8
    shell:
        """
        ./rules/scripts/run_vep_annotation.sh {wildcards.sample} {input.calls} {output.calls} {output.stats} {input.cache} {input.fasta} {input.gnomad_vcf} {input.dbsnp_vcf} {input.clinvar_vcf} {log} {params.cache_version} {params.species} {params.assembly}
        """

rule annotate_variants_svaba:
    input:
        calls="results/svaba/{sample}/{sample}_svaba.somatic.sv.vcf",
        cache="resources/vep/cache",
        fasta="resources/genome.fa",
        fai="resources/genome.fa.fai",
        gnomad_vcf="resources/gnomad.vcf.gz",
        dbsnp_vcf="resources/dbSNPv155_common.vcf.gz",
        clinvar_vcf="resources/clinvar.vcf.gz"
    output:
        calls="results/vep/{sample}/{sample}_svaba_somatic_sv_annotated.vcf",
        stats="results/vep/{sample}/{sample}_svaba_somatic_sv_annotated.html"
    params:
        species="homo_sapiens_gca009914755v4",
        assembly="T2T-CHM13v2.0",
        cache_version=106
    log:
        "logs/vep/{sample}/svaba_somatic_sv_annotated.log"
    threads: 8
    shell:
        """
        ./rules/scripts/run_vep_annotation.sh {wildcards.sample} {input.calls} {output.calls} {output.stats} {input.cache} {input.fasta} {input.gnomad_vcf} {input.dbsnp_vcf} {input.clinvar_vcf} {log} {params.cache_version} {params.species} {params.assembly}
        """

rule vcf_to_maf:
    input:
        vcf="results/vep/{sample}/{sample}_mutect_pass_variants_annotated.vcf",
        ref="resources/genome.fa"
    output:
        maf="results/vcf2maf/{sample}/{sample}_mutect_pass_variants_annotated.maf"
    params:
        vep_path="/staging/leuven/stg_00096/home/rdewin/system/miniconda/envs/VEP/bin",
        ncbi_build="T2T-CHM13v2.0",
    log:
        "logs/vcf2maf/{sample}_mutect_pass_variants.log"
    shell:
        """
        echo "Converting VCF to MAF..."
        perl {params.vep_path}/vcf2maf.pl \
            --input-vcf {input.vcf} \
            --output-maf {output.maf} \
            --tumor-id {wildcards.sample}_tumor \
            --normal-id {wildcards.sample}_normal \
            --vcf-tumor-id {wildcards.sample}_tumor \
            --vcf-normal-id {wildcards.sample}_normal \
            --ref-fasta {input.ref} \
            --ncbi-build {params.ncbi_build} \
            --inhibit-vep \
            &> {log}
        echo "VCF to MAF conversion completed."
        """

# rule annotate_with_annovar:
# 	input:
# 		vcf="path/to/input.vcf"
# 	output:
# 		annovar="path/to/output.annovar"
# 	params:
# 		db="path/to/annovar/database",
# 		protocol="gene,exac03,gnomad",  # Confirm protocols based on your needs
# 		operation="g,f,f",  # Confirm operations based on the protocols
# 		buildver="hs1"  # Confirm the build version
# 	shell:
# 		"""
# 		table_annovar.pl {input.vcf} {params.db} \
# 			-buildver {params.buildver} \
# 			-out {output.annovar} \
# 			-remove \
# 			-protocol {params.protocol} \
# 			-operation {params.operation} \
# 			-nastring . \
# 			-vcfinput
# 		"""

# rule annotate_with_vep:
# 	input:
# 		vcf="path/to/input.vcf"
# 	output:
# 		annotated_vcf="path/to/output.vep.vcf"
# 	params:
# 		vep_options="--cache --offline --everything",  # Adjust VEP options as needed
# 		species="homo_sapiens",  # Adjust the species as necessary
# 		assembly="GRCh38"  # Adjust the assembly version as necessary
# 	shell:
# 		"""
# 		vep --input_file {input.vcf} \
# 			--output_file {output.annotated_vcf} \
# 			--species {params.species} \
# 			--assembly {params.assembly} \
# 			{params.vep_options} \
# 			--vcf
# 		"""

# rule annotate_with_snpeff:
# 	input:
# 		vcf="path/to/input.vcf"
# 	output:
# 		annotated_vcf="path/to/output.snpeff.vcf"
# 	params:
# 		snpeff_genome="GRCh38.86"  # Adjust the genome version as necessary
# 	shell:
# 		"""
# 		snpEff {params.snpeff_genome} {input.vcf} > {output.annotated_vcf}
#     """
