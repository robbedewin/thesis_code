## Check allelically imbalanced and downregulated genes for nonsense mutations

airesults <- read.delim(file = "/staging/leuven/stg_00096/home/rdewin/ASE/results/allelic_imbalance_pooledsamples_vst.txt", sep = "\t", as.is = T)

airecurrences <- table(airesults[airesults$log2fc > 0, "gene_name"])
airecurrences <- sort(airecurrences[airecurrences > 1], decreasing = T)
airecurrences

aioccurrencespergene <- do.call(rbind, by(data = airesults, INDICES = airesults$gene_name, FUN = function(x) data.frame(samples = paste0(x$sample, collapse = ","),
                                                                                                                n_up = sum(x$log2fc > 0), n_down = sum(x$log2fc < 0), 
                                                                                                                updown = paste0(ifelse(x$log2fc > 0, "+", "-"), collapse = ","),
                                                                                                                gene_name = x$gene_name[1])))

write.table(x = aioccurrencespergene[, c("gene_name", "samples", "updown", "n_up", "n_down")],
            file = "/staging/leuven/stg_00096/home/rdewin/ASE/results/alloccurences_vst.txt", quote = F, sep = "\t", row.names = F)

