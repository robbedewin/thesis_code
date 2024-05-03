rule ascat:
    input:
        tumor_bam="results/recal/{sample}_dna_tumor_recal.bam",
        normal_bam="results/recal/{sample}_dna_normal_recal.bam"
    output:
        ascat_objects="results/ascat/{sample}_ASCAT_objects.Rdata",
        before_correction_tumour="results/ascat/plots/{sample}_Before_correction_{sample}_tumor.tumour.png",
        before_correction_germline="results/ascat/plots/{sample}_Before_correction_{sample}_tumor.germline.png",
        after_correction_tumour="results/ascat/plots/{sample}_After_correction_{sample}_tumor.tumour.png",
        after_correction_germline="results/ascat/plots/{sample}_After_correction_{sample}_tumor.germline.png",
        aspcf="results/ascat/plots/{sample}_{sample}_tumor.ASPCF.png"
    log:
        "logs/ascat/{sample}.log"
    threads: 8
    script:
        "scripts/ascat.R"


    