rule svaba:
    input:
        tumor_bam = "results/recal/{sample}/{sample}_dna_tumor_recal.bam",
        normal_bam = "results/recal/{sample}/{sample}_dna_normal_recal.bam",
    output:
        sv_somatic = "results/svaba/{sample}/{sample}_svaba.somatic.sv.vcf",
        sv_germline = "results/svaba/{sample}/{sample}_svaba.germline.sv.vcf",
        indel_somatic = "results/svaba/{sample}/{sample}_svaba.somatic.indel.vcf",
        indel_germline = "results/svaba/{sample}/{sample}_svaba.germline.indel.vcf",
        unfiltered_sv_somatic = "results/svaba/{sample}/{sample}_svaba.unfiltered.somatic.sv.vcf",
        unfiltered_sv_germline = "results/svaba/{sample}/{sample}_svaba.unfiltered.germline.sv.vcf",
        unfiltered_indel_somatic = "results/svaba/{sample}/{sample}_svaba.unfiltered.somatic.indel.vcf",
        unfiltered_indel_germline = "results/svaba/{sample}/{sample}_svaba.unfiltered.germline.indel.vcf",
        alignments_txt = "results/svaba/{sample}/{sample}_alignments.txt.gz",
        bps_txt = "results/svaba/{sample}/{sample}_bps.txt.gz",
        contigs_bam = "results/svaba/{sample}/{sample}_contigs.bam",
        discordant_txt = "results/svaba/{sample}/{sample}_discordant.txt.gz",
    params:
        ref = "resources/genome.fa",
        dbsnp = "resources/dbSNPv155_common.vcf.gz",
        #target = "chr17", -k {params.target} \
        blacklist = "resources/T2T.excluderanges.bed",
    log:
        log="logs/svaba/{sample}_svaba.log",
        full_log="logs/svaba/{sample}_svaba_full.log",
    threads: 8
    shell:
        """
        svaba run \
            -t {input.tumor_bam} -n {input.normal_bam} \
            -G {params.ref} \
            -D {params.dbsnp}\
            -B {params.blacklist} \
            -p {threads} \
            -a {wildcards.sample}_svaba \
            &> {log.log}
        
        # Move output files to specified locations
        mv {wildcards.sample}_svaba.log {log.full_log}

        mv {wildcards.sample}_svaba.svaba.somatic.sv.vcf {output.sv_somatic}
        mv {wildcards.sample}_svaba.svaba.germline.sv.vcf {output.sv_germline}
        mv {wildcards.sample}_svaba.svaba.somatic.indel.vcf {output.indel_somatic}
        mv {wildcards.sample}_svaba.svaba.germline.indel.vcf {output.indel_germline}
        mv {wildcards.sample}_svaba.svaba.unfiltered.somatic.sv.vcf {output.unfiltered_sv_somatic}
        mv {wildcards.sample}_svaba.svaba.unfiltered.germline.sv.vcf {output.unfiltered_sv_germline}
        mv {wildcards.sample}_svaba.svaba.unfiltered.somatic.indel.vcf {output.unfiltered_indel_somatic}
        mv {wildcards.sample}_svaba.svaba.unfiltered.germline.indel.vcf {output.unfiltered_indel_germline}

        mv {wildcards.sample}_svaba.alignments.txt.gz {output.alignments_txt}
        mv {wildcards.sample}_svaba.bps.txt.gz {output.bps_txt}
        mv {wildcards.sample}_svaba.contigs.bam {output.contigs_bam}
        mv {wildcards.sample}_svaba.discordant.txt.gz {output.discordant_txt}
        """