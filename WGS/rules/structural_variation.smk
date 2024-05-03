# rule delly_vcfgz:
#     input:
#         reference="resources/genome.fa",
#         alns=["mapped/a.bam"],
#         # optional
#         exclude="human.hg19.excl.tsv",
#     output:
#         "sv/calls.vcf.gz",
#     params:
#         extra="",  # optional parameters for delly (except -g, -x)
#     log:
#         "logs/delly.log",
#     threads: 2  # It is best to use as many threads as samples
#     wrapper:
#         "v3.7.0/bio/delly"


# Call structural variants with Delly
# Delly is a structural variant caller that uses paired-end mapping information to call deletions, duplications, inversions, and translocations.
rule delly_call:
    input:
        tumor="results/recal/{sample}/{sample}_dna_tumor_recal.bam",
        normal="results/recal/{sample}/{sample}_dna_normal_recal.bam",
        reference="resources/genome.fa",
        blacklist="resources/T2T.excluderanges.bed",
    output:
        "results/delly/call/{sample}/{sample}.bcf"
    log:
        "logs/delly/call/{sample}/{sample}.log"
    shell:
        """
        delly call \
            -x {input.blacklist} \
            -g {input.reference} \
            -o {output} \
            {input.tumor} \
            {input.normal} \
            &> {log} 
        """



rule delly_prefilter:
    input:
        bcf="results/delly/call/{sample}/{sample}.bcf",
        samples="config/samples.tsv"
    output:
        "results/delly/prefiltered/{sample}/{sample}.pre.bcf"
    log:
        "logs/delly/prefiltered/{sample}/{sample}.log"
    shell:
        """
        delly filter \
            -f somatic \
            -o {output} \
            -s {input.samples} \
            {input.bcf}
        """


rule delly_genotype:
    input:
        reference="resources/genome.fa",
        exclude="resources/T2T.excluderanges.bed",
        pre_bcf=expand("results/delly/prefiltered/{sample}/{sample}.pre.bcf", sample=unique_samples),
        bams=expand("results/recal/{sample}/{sample}_dna_{alias}_recal.bam", sample=unique_samples, alias=["tumor", "normal"])
    output:
        "results/delly/all.genotyped.bcf"
    shell:
        """
        delly call \
            -g {input.ref} \
            -v {input.pre_bcf} \
            -o {output} \
            -x {input.exclude} \
            {input.bams}
        """