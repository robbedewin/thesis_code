###!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### precancer ASCAT preprocessing
TUMOURIDX <- as.integer(args[1])
# TUMOURNAME <- "S094"


## analysis pipeline
# library(rslurm)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.GRCh38)
bsgnom <- BSgenome.Hsapiens.NCBI.GRCh38
# library(biomaRt)
library(ggplot2)
library(VGAM)
library(readr)
library(parallel)

source(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/1000Genomes_getAllelecounts.R")
source(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/utils.R")
source(file = "/staging/leuven/stg_00096/home/rdewin/WGS/rules/scripts/ase_nomatch_hg38_functions_refractoryALL.R")

## statics
REFANNOT <- "/camp/lab/vanloop/working/camp_pipeline_files/human/references/STAR-Fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021/ref_annot.gtf"
ALLELECOUNTER <- "/camp/lab/vanloop/working/demeulj/software/alleleCount-4.2.0/bin/alleleCounter"
ALLELESDIR <- "/camp/project/proj-vanloo/reference_files/human/references/Battenberg_hg38/1000G_loci_hg38/"
# JAVA <- "/srv/sw/eb/software/Java/1.8.0_221/bin/java"

GENOMEDIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2021_OConnor_refractory_T-ALL/data/WGS/mapped/"
RNADIR <- "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2021_OConnor_refractory_T-ALL/data/RNA-Seq/mapped/"
RNAREFGENOME <- "/camp/lab/vanloop/working/camp_pipeline_files/human/references/alignment/hs38DH/hs38DH.fa"

minMapQ <- 35
minBaseQ <- 20

# for cluster
NCORES <- 1

chrs <- c(1:22, "X")
# sampledf2 <- read.delim(file = "/camp/project/proj-emedlab-vanloo/jdemeul/projects/2019_fitzgibbon_ppAML/data/PoorRiskSamples.rna.and.dna.ondisk.tsv", as.is = T)
sampledf_full <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2021_OConnor_refractory_T-ALL/data/20210826_full_cohort.txt", as.is = T)
sampledf <- sampledf_full[, c("WGS_Pres_FASTQ", "WGS_Pres_FASTQ", "R_Pres_FASTQ")]
sampledf_rem <- sampledf_full[sampledf_full$WGS_Relapse_FASTQ != "", c("WGS_Relapse_FASTQ", "WGS_Relapse_FASTQ", "R_Relapse_FASTQ")]
colnames(sampledf) <- c("patient_id", "wgs_bam", "rna_bam")
colnames(sampledf_rem) <- c("patient_id", "wgs_bam", "rna_bam")

sampledf <- rbind(sampledf, sampledf_rem)
sampledf <- sampledf[sampledf$rna_bam != "", ]

sampledf$wgs_bam <- paste0(sampledf$wgs_bam, ".aln.recal.final.bam")
# sampledf$rna_bam <- paste0(sampledf$rna_bam, ".bam")


### read reference annotation
gc31 <- import.gff(con = file(REFANNOT), feature.type = "gene", colnames = c("gene_id", "gene_name", "gene_type"))


get_ase_ppaml <- function(sample_id, sample_id_rna) {
  
  # sample_id = sampledf$tumor_vial_id[1]
  # sample_id_rna <- sampledf[sampledf$patient_id %in% sample_id, "rna_bam"]

  TBAMFILE <- file.path(GENOMEDIR, sample_id, paste0(sample_id, ".aln.recal.final.bam"))
  TRNABAMFILE <- file.path(RNADIR, sample_id_rna, paste0(sample_id_rna, ".bam"))

  TOUTDIR <- paste0("/camp/project/proj-emedlab-vanloo/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/", sample_id)
  
  dir.create(path = TOUTDIR)
  
  if (!file.exists(file.path(TOUTDIR, paste0(sample_id, "_hetSNPs_nomatch.vcf.bgz")))) {
  ## get allelecounts for tumour exome
  mapply(FUN = alleleCount,
           locifile = paste0(ALLELESDIR, "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr", chrs, ".txt"),
           outfile = file.path(TOUTDIR, paste0(sample_id, "_alleleCounts_chr", chrs, ".txt")),
           MoreArgs = list(bam = TBAMFILE, min_baq = minBaseQ, min_maq = minMapQ))#,
           # mc.cores = NCORES, mc.preschedule = F)
  
  lapply(X = chrs, FUN = get_alleles_chr_nomatch, allelesdir = ALLELESDIR, countsdir = TOUTDIR, sample_id = sample_id)
         # , mc.cores = NCORES, mc.preschedule = T)
  combine_loci_nomatch(countsdir = TOUTDIR, sample_id = sample_id, bsgnom = BSgenome.Hsapiens.NCBI.GRCh38)
  } 
  
  if (!file.exists(file.path(TOUTDIR, paste0(sample_id, "_asereadcounts_nomatch.rtable")))) {
  samtoolsrhcmd <- paste0("/camp/apps/eb/software/SAMtools/1.13-GCC-10.2.0/bin/samtools",
                               " addreplacerg",
                               " -@ ", NCORES,
                               " -r '@RG\tID:1\tLB:lib1\tSM:", sample_id_rna, "\tPL:ILLUMINA'",
                               " -o ", paste0(dirname(TRNABAMFILE), "/", sample_id_rna, ".rehead.bam"),
                               " ", TRNABAMFILE)

  samtoolidxscmd <- paste0("/camp/apps/eb/software/SAMtools/1.13-GCC-10.2.0/bin/samtools",
                           " index",
                           " -@ ", NCORES,
                           " ", paste0(dirname(TRNABAMFILE), "/", sample_id_rna, ".rehead.bam"))
  
  ## get allelecounts for tumour transcriptome
  asereadcmd <- ASEReadCount(gatk.exe = "/camp/lab/vanloop/working/demeulj/software/gatk-4.2.0.0/gatk",
                             hetSNPvcf = file.path(TOUTDIR, paste0(sample_id, "_hetSNPs_nomatch.vcf.bgz")),
                             bamfile = paste0(dirname(TRNABAMFILE), "/", sample_id_rna, ".rehead.bam"),
                             refgenome = RNAREFGENOME,
                             outfile = file.path(TOUTDIR, paste0(sample_id, "_asereadcounts_nomatch.rtable")),
                             minBaseQ = minBaseQ, minMapQ = minMapQ)
  
  # actually run the constructed commands above
  print(paste0("Running ", samtoolsrhcmd))
  system(command = samtoolsrhcmd, wait = T)
  # system(command = paste0("mv ", TRNABAMFILE, " ", TRNABAMFILE, ".old;",
  #                         " mv ", TRNABAMFILE, ".bai ", TRNABAMFILE, ".bai.old;",
  #                         " mv ", dirname(TRNABAMFILE), "/", sample_id_rna, ".rehead.bam ", TRNABAMFILE), wait = T)

  print(paste0("Running ", samtoolidxscmd))
  system(command = samtoolidxscmd, wait = T)
  
  print(paste0("Running ", asereadcmd))
  system(command = asereadcmd, wait = T)
  }
  
  
  # if (!file.exists(file.path(TOUTDIR, paste0(sample_id, "_ase_out.txt")))) {
    ## compute statistics
    asedf <- compute_pvals_nomatch(toutdir = TOUTDIR, tsample = sample_id, exclude_bad_snps = F)
    
    
    ## Some QC and plotting
    # specifically set bitmapType to cairo, as this seems to get messed up during rslurm submission to the nodes
    # options(bitmapType="cairo")
    # print(getOption(x = "bitmapType"))
    # assess filtering
    
    p2 <- ggplot(data = asedf, mapping = aes(x = pval, fill = filter <= 0.01)) + geom_histogram(binwidth = 0.01) + scale_y_log10()
    ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_filter.png")), plot = p2, dpi = 300, width = 10, height = 7)
    
    # sum(asedf$padj < .05)
    # ggd.qqplot(asedf[asedf$pval > 0, "pval"])
    # ggd.qqplot(p.adjust(asedf[asedf$pval > 0, "pval"], method = "fdr"))
    
    p3 <- ggqq(asedf[asedf$pval > 0, "pval"])
    ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_QQ.png")), plot = p3, dpi = 300, width = 10, height = 7)
    
    ## Gene annotation
    # asedf <- merge(x = asedf, y = ase_annotate(asedf = asedf[asedf$pval <= 0.01, ], bsgnom = BSgenome.Hsapiens.NCBI.GRCh38, annot = gc31), all.x = T)
    asedf <- ase_annotate(asedf = asedf, bsgnom = BSgenome.Hsapiens.NCBI.GRCh38, annot = gc31)
    asedf$contig <- factor(x = asedf$contig, levels = c(1:22, "X"))
    asedf <- asedf[order(asedf$contig, asedf$position), ]
    write_tsv(x = asedf, file = file.path(TOUTDIR, paste0(sample_id, "_ase_out.txt")))
    
    # final plot
    p4 <- plot_ase_manhattan(asedf = asedf)
    ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_manhattan.png")), plot = p4, dpi = 300, width = 20, height = 6)
  # }

  return(NULL)
}

# debug(alleleCount)
# debug(get_ase_ppaml)
# mclapply(X = sampledf$sample_id, FUN = get_ase_ppaml, mc.preschedule = F, mc.cores = 16)
sample_id <- sampledf[TUMOURIDX, "patient_id"]
sample_id_rna <- sampledf[TUMOURIDX, "rna_bam"]

print(paste0("Running tumour ", sample_id, " with RNA id ", sample_id_rna, "."))

get_ase_ppaml(sample_id = sample_id, sample_id_rna = sample_id_rna)
# amlasejob <- slurm_apply(f = get_ase_ppaml, params = sampledf[,"sample_id", drop = F], jobname = "ase_ppaml_restart", nodes = 9, cpus_per_node = 6, add_objects = ls(),
#                           pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(), submit = T)
# print_job_status(amlasejob)
# cancel_slurm(amlasejob)
