###!/usr/bin/env Rscript
### reattemt to normalize log2fc values

if (!requireNamespace("DESeq2", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("DESeq2")
}
## get log2-fold expression values of T-ALL samples
library(DESeq2)

#countsdir <- snakemake@input[["counts_dir"]]
#sampleFiles <- list.files(path = countsdir, pattern = "ReadsPerGene.out.tab", recursive = TRUE)
#sampleName <- snakemake@wildcards[["sample"]]
#sampleTable <- data.frame(sampleName = sampleName,
#                          fileName = sampleFiles,
#                          condition = "T-ALL")


#setwd(countsdir)

countsdir <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/star"

setwd(countsdir)

sampleFiles <- list.files(path = countsdir, pattern = "ReadsPerGene.out.tab", recursive = T)
sampleName <- sapply(X = strsplit(x = sampleFiles, split = "/", fixed = T), FUN = "[", 1)
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = "T-ALL")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = countsdir,
                                       design= ~ 1)
# ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >= 100, ]

dds <- estimateSizeFactors(ddsHTSeq)
dds_vst <- vst(dds, blind = T)
# head(assay(dds_vst), 3)
# library(vsn)
# meanSdPlot(assay(dds_vst))
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))

# hist(assay(dds_vst)[3,])

dds_vstmeans <- rowMeans(x = assay(dds_vst))
vst_fc <- assay(dds_vst) - dds_vstmeans
plot(dds_vstmeans, vst_fc[,1])
plot(log10(rowMeans(counts(dds, normalized=TRUE))+1), vst_fc[,1])

resdf <- as.data.frame(counts(dds, normalized=TRUE))

l2fcdf <- as.data.frame(vst_fc)

#gene_ids_names <- read.delim(file = "/camp/lab/vanloop/working/camp_pipeline_files/human/references/STAR-Fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021/ref_annot.gtf.gene_spans", as.is = T, header = F)
gene_ids_names <- read.delim(file = "/staging/leuven/stg_00096/references/CTAT_Resources/T2T-CHM13_CTAT_lib_Feb162023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans", as.is = T, header = F)
rownames(gene_ids_names) <- gene_ids_names$V1


l2fcdf$gene_name <- gene_ids_names[rownames(l2fcdf), "V6"]
l2fcdf$mean_expression <- 2^rowMeans(x = log2(resdf+1))

resdf$gene_name <- gene_ids_names[rownames(resdf), "V6"]

write.table(x = l2fcdf, file = "/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/RNAlog2fc_vst.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = as.data.frame(assay(dds_vst)), file = "/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/RNAcounts_vst.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = resdf, file = "/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/RNAcounts_normalised_T-ALL.txt", quote = F, sep = "\t", row.names = T, col.names = T)


if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    BiocManager::install("GenomicFeatures")
}

## combine p-values (of powered SNP loci) per gene
library(GenomicFeatures)
library(ggplot2)

### functions

combine_pvals <- function(ase_pergene) {
  # ase_pergene <- ase_results_annot[4:14, ]
  outdf <- data.frame(contig = ase_pergene[1, "contig"], positions = paste0(unique(ase_pergene$position), collapse = ","), pcombined = 1, gene = ase_pergene[1, "gene"], stringsAsFactors = F, power = T)
  
  is_duplicated <- duplicated(ase_pergene$position)
  has_power <- ase_pergene$filter <= 0.01
  ase_pergene <- ase_pergene[!is_duplicated & has_power, ]
  
  if (nrow(ase_pergene) > 1) {
    outdf$pcombined <- fishersMethod(x = ase_pergene$pval)
  } else if (nrow(ase_pergene) == 1) {
    outdf$pcombined <- ase_pergene$pval
  } else {
    outdf$power <- F
  }
  return(outdf)
}


fishersMethod <- function(x) {
  pchisq(q = -2 * sum(log(x)), df = 2*length(x), lower.tail = F)
}


plot_imbalance_expression <- function(imbalancedf) {
  imbalancedf <-  imbalancedf[order(imbalancedf$mean_expression, decreasing = F), ]
  labeldf <- data.frame(pos = unlist(lapply(X = 10^(0:4), FUN = function(x) sum(imbalancedf$mean_expression < x))), expr = 10^(0:4), stringsAsFactors = F)
  
  outdf_bak <- imbalancedf
  imbalancedf <- imbalancedf[!grepl(pattern = "^HLA.*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^IG[HLK].*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^TR[ABDG][VCDJ].*", x = imbalancedf$gene_name, perl = T), ]

  imbalancedf$notes <- ifelse(imbalancedf$padj > 0.05, "nonsig", 
                              ifelse(imbalancedf$log2fc >= 1, "up",
                                     ifelse(imbalancedf$log2fc <= -.73, "down", "nonsig")))
  
  p1 <- ggplot(data = imbalancedf, mapping = aes(x = 1:nrow(imbalancedf), y = -sign(log2fc)*log10(pcombined)))
  p1 <- p1 + geom_point(mapping = aes(colour = notes, size = abs(log2fc)), 
                        show.legend = F, alpha = .4)
  p1 <- p1 + geom_hline(yintercept = c(-1,1)*-log10(max(imbalancedf[imbalancedf$padj < .05, "pcombined"])), linetype = "dashed", colour = "grey") +
    geom_text(data = imbalancedf[imbalancedf$notes != "nonsig", ], mapping = aes(x = which(imbalancedf$notes != "nonsig"), y = -sign(log2fc)*log10(pcombined), label = gene_name), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(imbalancedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_y_continuous(breaks = seq(-10,10,2), oob = scales::squish, limits = c(-10,10))
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$expr, name = "mean expression (normalised)")
  # p1 <- p1 + scale_color_brewer(type = "div", palette = "RdBu", direction = -1)
  p1 <- p1 + scale_color_manual(values = c(nonsig = "#e0e0e0", up = "#ef8a62", down = "#67a9cf"))
  p1 <- p1 + scale_size_continuous(range = c(1,7.5))
  p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = -90)) + labs(x = NULL)
  return(p1)  
}
### end functions

# generate library of all Hs exons (should be the covered regions).
gtffile <- "/staging/leuven/stg_00096/references/CTAT_Resources/T2T-CHM13_CTAT_lib_Feb162023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
# seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb)))
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))
seqlevelsStyle(hsexondb) <- "Ensembl"

#sampledf <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2021_OConnor_refractory_T-ALL/data/20210826_full_cohort.txt", as.is = T)
#sampledf <- sampledf[sampledf$R_Pres_FASTQ != "", ]

# add in the log2-fold change data and actual gene names
l2fcfile <- "/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/RNAlog2fc_vst.txt"
l2fcdf <-  read.delim(file = l2fcfile, as.is = T)


for (SAMPLEID in sampleName) {
#for (i in 1:nrow(sampledf)) {
  # i <- 1
  #SAMPLEID <- sampledf[i, "R_Pres_FASTQ"]
  
  ## read a results file
  # if (sampledf[i, "cell_line"]) {
  #   ase_resultsfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/cell_lines/", SAMPLEID, "/", SAMPLEID, "_ase_out.txt")
  # } else {
    #WGSID <- sampledf[i, "WGS_Pres_FASTQ"]
    ase_resultsfile <- paste0("/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/", SAMPLEID, "/", SAMPLEID, "_ase_out.txt")
  # }
  outDir <- dirname(ase_resultsfile)

# Create the directory if it does not exist
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  ase_results <- read.delim(file = ase_resultsfile, as.is = T)
  
  
  # make results into GRanges object, identify all exonic SNPs and create new df with all of these (contains duplicate SNPs)
  asegr <- GRanges(seqnames = ase_results$contig, ranges = IRanges(start = ase_results$position, end = ase_results$position))
  annothits <- findOverlaps(query = asegr, subject = hsexondb)
  # in one case, there were two genes using the same exon ... this just takes the first
  hitgenes <- sapply(mcols(hsexondb[subjectHits(annothits)])$gene_id, FUN = function(x) x[[1]])
  ase_results_annot <- data.frame(ase_results[queryHits(annothits), colnames(ase_results) != "gene"], gene = hitgenes, stringsAsFactors = F)
  
  
  # create output dataframe with combined p-value per gene + adjust for multiple testing
  outdf <- do.call(rbind, by(data = ase_results_annot, INDICES = ase_results_annot$gene, FUN = combine_pvals))
  outdf$padj <- 1
  outdf[outdf$power, "padj"] <- p.adjust(p = outdf[outdf$power, "pcombined"], method = "fdr")
  # outdf$gene_name <- gene_ids_names[outdf$gene, "V2"]
  
  # add in the log2-fold change data and actual gene names
  outdf[, c("log2fc", "mean_expression", "gene_name")] <- l2fcdf[outdf$gene, c(grep(pattern = paste0(sub(pattern = "-", replacement = ".", SAMPLEID), "$"), x = colnames(l2fcdf), value = T), "mean_expression", "gene_name")]
  
  # format
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  # if (sampledf[i, "cell_line"]) {
  #   outfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/cell_lines/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.txt")
  # } else {
    outfile <- paste0("/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.txt")
  # }
  
  # plot
  p1 <- plot_imbalance_expression(imbalancedf = outdf)
  # if (sampledf[i, "cell_line"]) {
  #   plotfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/cell_lines/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.png")
  # } else {
    plotfile <- paste0("/staging/leuven/stg_00096/home/rdewin/RNA/results/ASE/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.png")
  # }
  ggsave(filename = plotfile, plot = p1, dpi = 300, width = 15, height = 6)
}
