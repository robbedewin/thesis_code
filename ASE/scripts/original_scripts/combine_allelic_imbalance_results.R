## check recurrence of allelically imbalanced + up/downregulated genes

resultfiles <- list.files(path = "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T)
ai_results <- lapply(X = resultfiles, FUN = read.delim, as.is = T)
# View(ai_results[[1]])

aidf <- do.call(rbind, ai_results)
aidf$sample <- rep(sub(pattern = "_imbalance_expression_vst.txt", replacement = "", x = basename(resultfiles)), sapply(X = ai_results, nrow))

aidf <- aidf[aidf$padj <= .05 & (aidf$log2fc >= 1 | aidf$log2fc <= -0.73), ]
aidf <- aidf[aidf$contig != "X" & !grepl(aidf$gene_name, pattern = "HLA-*") &
               !grepl(aidf$gene_name, pattern = "^IG[HLK].*", perl = T) &
               !grepl(aidf$gene_name, pattern = "^TR[ABDG][VCDJ].*", perl = T), ]

recurrent_genes <- sort(table(aidf$gene_name), decreasing = T)
head(recurrent_genes, n = 25)

write.table(x = aidf, file = "/camp/project/proj-vanloo/analyses/jdemeul/projects/2021_OConnor_refractory_T-ALL/results/ASE/20220320_allelic_imbalance_pooledsamples_vst.txt", sep = "\t", quote = F, row.names = F, col.names = T)
