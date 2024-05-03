## create final plot to show recurrences and outlying log2fc

# get sample IDs mappings
sampledf_full <- read.delim(file = "/camp/lab/vanloop/working/demeulj/projects/2021_OConnor_refractory_T-ALL/data/20210826_full_cohort.txt", as.is = T)
sampledf <- sampledf_full[sampledf_full$R_Pres_FASTQ != "", ]

# get log2fc matrix (VST)
l2fcfile <- "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/20220320_RNAlog2fc_vst.txt"
l2fcdf <-  read.delim(file = l2fcfile, as.is = T)

# reformat column names
# colnames(l2fcdf)[sapply(X = paste0(sub(pattern = "-", replacement = ".", sampledf$sampleid), "$"), FUN = grep, x = colnames(l2fcdf), simplify = T)] <- c(paste0("WES_", sampledf[!sampledf$cell_line, "t_wes_id"]), sub(pattern = "-", replacement = ".", x = sampledf[sampledf$cell_line, "sampleid"]))
# colnames(l2fcdf)[grep(pattern = "^R", x = colnames(l2fcdf))] <- sampledf[match(x = grep(pattern = "^R", x = colnames(l2fcdf), value = T), table = sampledf$R_Pres_FASTQ), "WGS_Pres_FASTQ"]


# get allelic imbalance pooled samples file
airesultsfile <- "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/20220320_alloccurences_vst.txt"
airesults <- read.delim(file = airesultsfile, sep = "\t", as.is = T)
airesults_sub <- airesults[airesults$n_up > 2, ]

# splitsamplestring <- strsplit(airesults_sub$samples, split = ",", fixed = T)
# imbalanced_gene_sample_combinations <- paste0(rep(x = airesults_sub$gene_name, lengths(splitsamplestring)), "_", unlist(splitsamplestring))
# 
l2fcdf_sub <- l2fcdf[l2fcdf$gene_name %in% airesults_sub$gene_name, ]
# l2fcdf_sub$exprrank <- (1:nrow(l2fcdf_sub))[order(l2fcdf_sub$mean_expression, decreasing = F)]
l2fcdf_sub$gene_name <- factor(x = l2fcdf_sub$gene_name, levels = l2fcdf_sub[order(l2fcdf_sub$mean_expression, decreasing = F), "gene_name"])
# subset log2fc mat to the represented genes (maybe those that occur at least twice to get overview on X)



### MOD, include ALL AI cases irrespective of log2fc
resultfiles <- c(list.files(path = "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T))

airesults_all <- lapply(X = resultfiles, FUN = read.delim, as.is = T)
# View(ai_results[[1]])

airesults_all_df <- do.call(rbind, airesults_all)
airesults_all_df$sample_WGS <- sub(pattern = "-", replacement = ".", rep(sub(pattern = "_imbalance_expression_vst.txt", replacement = "", x = basename(resultfiles)), sapply(X = airesults_all, nrow)))
airesults_all_df$sample <- sampledf[match(x = airesults_all_df$sample_WGS, table = sampledf$WGS_Pres_FASTQ), "R_Pres_FASTQ"]

airesults_all_df <- airesults_all_df[airesults_all_df$padj <= .05, ]
airesults_all_df <- airesults_all_df[airesults_all_df$contig != "X" & !grepl(airesults_all_df$gene_name, pattern = "HLA-*") &
                                       !grepl(airesults_all_df$gene_name, pattern = "^IG[HLK].*", perl = T) &
                                       !grepl(airesults_all_df$gene_name, pattern = "^TR[ABDG][VCDJ].*", perl = T), ]

imbalanced_gene_sample_combinations <- paste0(airesults_all_df$gene_name, "_", airesults_all_df$sample)
###



library(reshape2)
l2fcdf_melt <- melt(data = l2fcdf_sub, id.vars = c("gene_name", "mean_expression"), variable.name = "sample_id", value.name = "l2fc")
l2fcdf_melt$is_ai <- paste0(l2fcdf_melt$gene_name, "_", l2fcdf_melt$sample_id) %in% imbalanced_gene_sample_combinations

# plot mean counts (rank axis) vs log2fc_vst for all samples for those genes
# simple dot for each, larger dot for samples in panel and fill for imbalanced ones
library(ggplot2)

p1 <- ggplot(data = l2fcdf_melt, mapping = aes(x = gene_name, y = l2fc))
# p1 <- p1 + geom_point(alpha = .5, shape = 16, size = 1.5, stroke = 0)
p1 <- p1 + geom_violin(scale = "width", alpha =.5)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$sample_id %in% c("R005", "R007", "R113"), ], colour = "green", alpha = .5, size = 2, stroke = 0)
p1 <- p1 + geom_point(data = l2fcdf_melt[!l2fcdf_melt$sample_id %in% c("R005", "R007", "R113"), ], alpha = .5, size = 2, stroke = 0)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$is_ai, ], colour = "red", alpha = .6, shape = 1, size = 2, stroke = 1)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
# p1
ggsave(filename = "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/20220320_l2fc_vst_AIrecurrenceGreaterThan2.png", plot = p1, dpi = 300, width = 30, height = 6)
write.table(x = l2fcdf_melt, file = "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/20220320_l2fc_vst_AIrecurrenceGreaterThan2_table.txt", sep = "\t", col.names = T, row.names = F, quote = F)
