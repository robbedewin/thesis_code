rule modify_bam_header:
    input:
        bam="results/recal/{sample}/{sample}_dna_{alias}_recal.bam"
    output:
        temp_bam=temp("results/recal/{sample}/{sample}_dna_{alias}_recal_temp.bam")
    shell:
        """
        samtools view -H {input.bam} > header.txt
        sed -i 's/SM:{wildcards.sample}/SM:{wildcards.sample}_{wildcards.alias}/g' header.txt
        samtools reheader header.txt {input.bam} > {output.temp_bam}
        mv {output.temp_bam} {input.bam}
        """