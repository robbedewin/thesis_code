rule liftover_chm13_to_grch37_bam:
    input:
        original="results/recal/{sample}/{sample}_dna_{alias}_recal.bam",
        chain="resources/liftover/chm13v2-to-grch37.chain"
    output:
        lifted="results/liftover/{sample}/{sample}_dna_{alias}_recal.lifted_grch37.bam",
    shell:
        "CrossMap.py bam {input.chain} {input.original} {output.lifted}"

rule liftover_chm13_to_grch37_vcf:
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