rule haplotype_caller:
    input:
        recal_bam="results/recal/{sample}_dna_{alias}_recal.bam",
        genome="resources/genome.fa"
    output:
        gvcf=temp("results/gvcf/{sample}_dna_{alias}.g.vcf"),
        gvcf_idx=temp("results/gvcf/{sample}_dna_{alias}.g.vcf.idx")
    log:
        haplo_log="logs/gatk/{sample}_dna_{alias}_haplotypecaller.log"
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
            -R {input.genome} \
            -I {input.recal_bam} \
            -O {output.gvcf} \
            -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
            -ERC GVCF &> {log.haplo_log}
        """

rule combine_gvcfs:
    input:
        gvcf_list="resources/gvcf_list.txt", # A text file listing all GVCF files to be combined
        genome="resources/genome.fa"
    output:
        combined_gvcf="results/gvcf/combined.g.vcf",
        combined_gvcf_idx="results/gvcf/combined.g.vcf.idx"
    log:
        combine_log="logs/gatk/combine_gvcfs.log"
    shell:
        """
        gatk --java-options "-Xmx4g" CombineGVCFs \
            -R {input.genome} \
            --variant {input.gvcf_list} \
            -O {output.combined_gvcf} &> {log.combine_log}
        """

rule genotype_gvcfs:
    input:
        combined_gvcf="results/gvcf/combined.g.vcf",
        genome="resources/genome.fa"
    output:
        genotyped_vcf="results/vcf/{sample}_combined_genotyped.vcf",
        genotyped_vcf_idx="results/vcf/{sample}_combined_genotyped.vcf.idx"
    log:
        genotype_log="logs/gatk/genotype_gvcfs.log"
    shell:
        """
        gatk --java-options "-Xmx4g" GenotypeGVCFs \
            -R {input.genome} \
            -V {input.combined_gvcf} \
            -O {output.genotyped_vcf} \
            -G StandardAnnotation -G AS_StandardAnnotation &> {log.genotype_log}
        """
rule variant_recalibrator_snps:
    input:
        vcf="results/vcf/{sample}_combined_genotyped.vcf",
        ref="resources/genome.fa",
        hapmap="resources/hapmap.vcf",
        omni="resources/omni.vcf",
        thousandG="resources/1000G.vcf",
        dbsnp="resources/dbsnp.vcf"
    output:
        recal="results/vcf/{sample}_snps.recal",
        tranches="results/vcf/{sample}_snps.tranches",
        plots="results/vcf/{sample}_snps_plots.R"
    log:
        vr_log_snps="logs/gatk/{sample}_variant_recalibrator_snps.log"
    shell:
        """
        gatk VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
            -mode SNP \
            --tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} \
            &> {log.vr_log_snps}
        """
rule variant_recalibrator_indels:
    input:
        vcf="results/vcf/{sample}_combined_genotyped.vcf",
        ref="resources/genome.fa",
        mills="resources/Mills_and_1000G_gold_standard.indels.vcf",
        dbsnp="resources/dbsnp.vcf"
    output:
        recal="results/vcf/{sample}_indels.recal",
        tranches="results/vcf/{sample}_indels.tranches",
        plots="results/vcf/{sample}_indels_plots.R"
    log:
        vr_log_indels="logs/gatk/{sample}_variant_recalibrator_indels.log"
    shell:
        """
        gatk VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
            -mode INDEL \
            --tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} \
            &> {log.vr_log_indels}
        """
rule apply_vqsr_snps:
    input:
        vcf="results/vcf/{sample}_combined_genotyped.vcf",
        recal="results/vcf/{sample}_snps.recal",
        tranches="results/vcf/{sample}_snps.tranches",
        ref="resources/genome.fa"
    output:
        vcf="results/vcf/{sample}_recalibrated_snps.vcf"
    log:
        apply_log_snps="logs/gatk/{sample}_apply_vqsr_snps.log"
    shell:
        """
        gatk ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranches} \
            -mode SNP \
            --truth-sensitivity-filter-level 99.0 \
            -O {output.vcf} \
            &> {log.apply_log_snps}
        """
rule apply_vqsr_indels:
    input:
        vcf="results/vcf/{sample}_recalibrated_snps.vcf",
        recal="results/vcf/{sample}_indels.recal",
        tranches="results/vcf/{sample}_indels.tranches",
        ref="resources/genome.fa"
    output:
        vcf="results/vcf/{sample}_recalibrated_snps_indels.vcf"
    log:
        apply_log_indels="logs/gatk/{sample}_apply_vqsr_indels.log"
    shell:
        """
        gatk ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranches} \
            -mode INDEL \
            --truth-sensitivity-filter-level 99.0 \
            -O {output.vcf} \
            &> {log.apply_log_indels}
        """
