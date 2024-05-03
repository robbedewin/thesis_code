rule link_files:
    output:
        ref_genome="resources/genome.fa",
        ref_genome_index="resources/genome.fa.fai",
        ref_genome_dict="resources/genome.dict",
        annotation="resources/annotation.gtf",
        bwa_index1="resources/genome.fa.0123",
        bwa_index2="resources/genome.fa.amb",
        bwa_index3="resources/genome.fa.ann",
        bwa_index4="resources/genome.fa.bwt.2bit.64",
        bwa_index5="resources/genome.fa.pac",

    shell:
        """
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa {output.ref_genome}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.fai {output.ref_genome_index}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.dict {output.ref_genome_dict}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/annotation.gtf {output.annotation}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.0123 {output.bwa_index1}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.amb {output.bwa_index2}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.ann {output.bwa_index3}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.bwt.2bit.64 {output.bwa_index4}
        ln -s /staging/leuven/stg_00096/home/rdewin/WGS/resources/genome.fa.pac {output.bwa_index5}
        """

rule minimap2_index:
    input:
        ref_genome="resources/genome.fa"
    output:
        index="resources/genome.mmi"
    shell:
        """
        minimap2 -d {output.index} {input.ref_genome}
        """