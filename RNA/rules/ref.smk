rule link_files:
    output:
        ref_genome="resources/genome.fa",
        ref_genome_index="resources/genome.fa.fai",
        ref_genome_dict="resources/genome.dict",
        annotation="resources/annotation.gtf",
        ctat_lib=directory("resources/CTAT_lib.plug-n-play"),
        ctat_lib_source=directory("resources/CTAT_lib.source"),
        
    shell:
        """
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa {output.ref_genome}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.fai {output.ref_genome_index}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.dict {output.ref_genome_dict}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf {output.annotation}
        ln -s /staging/leuven/stg_00096/references/CTAT_Resources/T2T-CHM13_CTAT_lib_Feb162023.plug-n-play {output.ctat_lib}
        ln -s /staging/leuven/stg_00096/references/CTAT_Resources/T2T-CHM13_CTAT_lib_Feb162023.source {output.ctat_lib_source}
        """

# Creates an index of the reference genome for the STAR RNA-seq aligner
rule STAR_index:
    input:
        genome="resources/genome.fa",
        annotation="resources/annotation.gtf"
    output:
        directory("resources/star_genome")
    threads: 4
    log:
        "logs/reference/star_index.log"
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation} \
            --sjdbOverhang 100
        && mv {output}/Log.out logs/reference/star_index.log
        """

