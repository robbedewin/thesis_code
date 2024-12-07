## Create final plot to show recurrences and outlying log2fc

# Get log2fc matrix (VST)
l2fcfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/expression_data/RNAlog2fc_vst.txt"
l2fcdf <- read.delim(file = l2fcfile, as.is = T)


## Get sample IDs 
matchedSamples <- c("P011", "P013", "P016", "P017", "P018", "P019", "P020", "P022", "P023", "P024", "P026", "P028", "P029", "P033", "P037", "P041", "P057", "P058", "P059", "P060", "P061", "P064", "P065", "P066", "P086", "P103", "P105")

## Get allelic imbalance pooled samples file
airesultsfile <- "/staging/leuven/stg_00096/home/rdewin/ASE/results/alloccurences_vst.txt"
airesults <- read.delim(file = airesultsfile, sep = "\t", as.is = T)
airesults_sub <- airesults[airesults$n_up > 2, ]


l2fcdf_sub <- l2fcdf[l2fcdf$gene_name %in% airesults_sub$gene_name, ]

l2fcdf_sub$gene_name <- factor(x = l2fcdf_sub$gene_name, levels = l2fcdf_sub[order(l2fcdf_sub$mean_expression, decreasing = F), "gene_name"])


library(reshape2)

imbalanced_gene_sample_combinations <- paste0(airesults_sub$gene_name, "_", airesults_sub$samples)

l2fcdf_melt <- melt(data = l2fcdf_sub, id.vars = c("gene_name", "mean_expression"), variable.name = "sample_id", value.name = "l2fc")
l2fcdf_melt$is_ai <- paste0(l2fcdf_melt$gene_name, "_", l2fcdf_melt$sample_id) %in% imbalanced_gene_sample_combinations

# plot mean counts (rank axis) vs log2fc_vst for all samples for those genes
# simple dot for each, larger dot for samples in panel and fill for imbalanced ones
library(ggplot2)

p1 <- ggplot(data = l2fcdf_melt, mapping = aes(x = gene_name, y = l2fc))
p1 <- p1 + geom_point(alpha = .5, shape = 16, size = 1.5, stroke = 0)
p1 <- p1 + geom_violin(scale = "width", fill = "grey", alpha = .5)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$is_ai, ], aes(x = gene_name, y = l2fc, fill = is_ai), shape = 21, size = 2, stroke = 0)
p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggsave(filename = "/staging/leuven/stg_00096/home/rdewin/ASE/results/l2fc_vst_AIrecurrenceGreaterThan2.png", plot = p1, width = 10, height = 5, dpi = 300)
write.table(x = l2fcdf_sub, file = "/staging/leuven/stg_00096/home/rdewin/ASE/results/l2fc_vst_AIrecurrenceGreaterThan2.txt", quote = F, sep = "\t", row.names = F)
