# Downloads and decompresses the reference genome
rule get_reference:
    output:
        "resources/genome.fa"
    log:
        "logs/reference/reference.log"
    shell:
        """
        wget -O {output}.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz
        gunzip {output}.gz
        """

#Generates FASTA index of the reference genome using samtools
rule genome_faidx:
    input:
        "resources/genome.fa"
    output:
        "resources/genome.fa.fai"
    log:
        "logs/reference/faidx.log"
    shell:
        """
        samtools faidx {input} -o {output} 2> {log}
        """

#Creates a dictionary for the reference genome, required by certain analysis tools
rule genome_dict:
    input:
        "resources/genome.fa"
    output:
        "resources/genome.dict"
    log:
        "logs/reference/dict.log"
    shell:
        """
        samtools dict {input} > {output} 2> {log}
        """

# Downloads gene annotations and converts GFF3 to GTF format
rule get_annotation:
    output:
        gtf="resources/annotation.gtf",
        gff3="resources/chm13v2.0.gff3"
    log:
        "logs/reference/annotation.log"
    shell:
        """
        wget -O {output.gff3}.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz
        gunzip {output.gff3}.gz
        gffread {output.gff3} -T -o {output.gtf}
        """

# Filters common variants from dbSNP for the reference genome
rule get_common_dbSNP:
    output:
        vcf="resources/dbSNPv155.vcf.gz",
        vcf_index="resources/dbSNPv155.vcf.gz.tbi",
        vcf_common="resources/dbSNPv155_common.vcf.gz",
        vcf_common_index="resources/dbSNPv155_common.vcf.gz.tbi"
    log:
        "logs/reference/get_dbSNP.log",
    shell:
        """
        (
        wget -O {output.vcf} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_dbSNPv155.vcf.gz
        bcftools index -t {output.vcf} -o {output.vcf_index}
        bcftools view -i INFO/COMMON=1 {output.vcf} -O z -o {output.vcf_common}
        bcftools index -t {output.vcf_common} -o {output.vcf_common_index}
        ) >> {log} 2>&1
        """

# Indexes the genome for BWA-MEM2, facilitating fast and accurate alignments
rule bwa_index:
    input:
        "resources/genome.fa"
    output:
        multiext("resources/genome.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        "logs/reference/bwa_index.log"
    shell:
        """
        bwa-mem2 index {input} 2> {log}
        """

# Downloads germline variant data from gnomAD for variant analysis
rule germline_resource:
    output:
        vcf="resources/gnomad.vcf.gz",
        vcf_tbi="resources/gnomad.vcf.gz.tbi",
    log:
        "logs/reference/germline_resource.log"
    shell:
        """
        wget -O {output.vcf} https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz
        wget -O {output.vcf_tbi}  https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz.tbi
        """

rule germline_resource_chr:
    input:
        vcf="resources/gnomad.vcf.gz",
    output:
        vcf_chr="resources/gnomad_chr.vcf.gz",
        #vcf_chr_tbi="resources/gnomad_chr.vcf.gz.tbi",
    log:
        "logs/reference/germline_resource_chr.log",
    shell:
        """
        # Extract the header
        bcftools view -h {input.vcf} > header.txt

        # Modify the header to rename the contigs and write the modified header to a new file
        awk -F'=' '/^##contig=<ID=/{{sub(/ID=/, "ID=chr", $0)}}1' header.txt > modified_header.txt
        
        # Replace the header of the VCF file with the modified header and compress the output
        bcftools reheader -h modified_header.txt -o temp.vcf {input.vcf}
        bgzip -c -f temp.vcf > {output.vcf_chr}

        # Generate a new .tbi file
        # tabix -p vcf {output.vcf_chr}
        """

rule filter_gnomad_vcf:
    input:
        vcf_chr="resources/gnomad_chr.vcf.gz"
    output:
        filtered_vcf="resources/filtered_gnomad.vcf.gz"
    shell:
        """
        bcftools view -i 'INFO/AF>0.01 & INFO/AF<=0.2' -m2 -M2 -v snps {input.vcf_chr} -Oz -o {output.filtered_vcf}
        bcftools index -t {output.filtered_vcf}
        """

rule create_interval_list:
    input:
        filtered_vcf="resources/filtered_gnomad.vcf.gz"
    output:
        interval_list="resources/filtered_gnomad.interval_list"
    shell:
        """
        bcftools query -f'%CHROM\t%POS0\t%POS\t%ID\n' {input.filtered_vcf} > {output.interval_list}
        """

# Downloads the T2T-CHM13 exclusion regions BED file
rule get_exclude_regions:
    output:
        exclusion_bed="resources/T2T.excluderanges.bed"
    shell:
        """
        gdown --id 1p8qGWHzv8ayud1g_Vx9QDYhprlbWWjKD -O {output.exclusion_bed}
        """
#Downloaded from https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT 


# Downloads the ClinVar database for variant interpretation
rule get_clinvar:
    output:
        vcf="resources/clinvar.vcf.gz",
    log:
        "logs/reference/get_clinvar.log"
    shell:
        """
        wget -O {output.vcf} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz
        tabix -p vcf {output.vcf}
        &> {log}
        """


# Downloads the T2T-CHM13 chain files for liftover from GRCh38 and GRCh37
rule get_chain:
    output:
        chain_grch38_chm13="resources/liftover/grch38-to-chm13v2.chain",
        chain_chm13_grch38="resources/liftover/chm13v2-to-grch38.chain",
        chain_grch37_chm13="resources/liftover/grch37-to-chm13v2.chain",
        chain_chm13_grch37="resources/liftover/chm13v2-to-grch37.chain",
    shell:
        """
        wget -O {output.chain_grch38_chm13} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
        wget -O {output.chain_chm13_grch38} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
        wget -O {output.chain_grch37_chm13} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/hg19-chm13v2.chain
        wget -O {output.chain_chm13_grch37} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-hg19.chain
        """

# Liftover GRCh38 variants to CHM13v2
rule liftover:
    input:
        original="resources/liftover/grch38.vcf",
        chain="resources/grch38-to-chm13v2.chain"
    output:
        lifted="resources/liftover/chm13.lifted.vcf",
        unlifted="resources/liftover/chm13.unlifted.vcf"
    params:
        genome="resources/genome.fa" 
    shell:
        """
        picard LiftoverVcf \
            I={input.original} \
            O={output.lifted} \
            CHAIN={input.chain} \
            REJECT={output.unlifted} \
            R={params.genome}
        """

rule gridss_pon:
    output:
        "resources/gridss/pon_gridds.tar.gz"
    shell:
        """
        wget -O {output} "https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources%2FDNA-Resources&files=hmf_pipeline_resources.38_v5.31.gz"
        tar -xvzf {output} -C resources/gridss
        """

#Get the download link from the website https://www.openbioinformatics.org/annovar/annovar_download_form.php
rule get_annovar:
    output:
        directory("resources/annovar")
    shell:
        """
        wget -O annovar.tar.gz http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
        mkdir -p {output}
        tar -xvzf annovar.tar.gz -C {output} --strip-components 1
        rm annovar.tar.gz
        """

rule download_databases:
    input:
        annovar_dir="resources/annovar"
    output:
        touch("resources/annovar_humandb_download_complete")
    shell:
        """
        mkdir -p {input.annovar_dir}/humandb
        #{input.annovar_dir}/annotate_variation.pl -buildver hs1 -downdb -webfrom annovar refGene {input.annovar_dir}/humandb || true
        #{input.annovar_dir}/annotate_variation.pl -buildver hs1 -downdb -webfrom annovar cytoBand {input.annovar_dir}/humandb
        {input.annovar_dir}/annotate_variation.pl -buildver hs1 -downdb -webfrom annovar gnomad_genome {input.annovar_dir}/humandb
        
        # Additional steps for generating FASTA files
        {input.annovar_dir}/annotate_variation.pl --buildver hs1 --downdb seq {input.annovar_dir}/humandb/hs1_seq
        {input.annovar_dir}/retrieve_seq_from_fasta.pl {input.annovar_dir}/humandb/hs1_refGene.txt -seqdir {input.annovar_dir}/humandb/hs1_seq -format refGene -outfile {input.annovar_dir}/humandb/hs1_refGeneMrna.fa
        
        touch {output}
        """