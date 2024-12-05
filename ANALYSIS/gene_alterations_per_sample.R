# Extract gene names and sample IDs from MAF object
maf_genes_samples <- maf_combined %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  dplyr::group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()


# Extract gene names and sample IDs from copy number GRanges object
copy_number_genes_samples <- as.data.frame(cnv_annotated) %>%
  dplyr::select(gene_name, Sample) %>%
  dplyr::group_by(gene_name, Sample) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()


# Extract gene names and sample IDs from structural variants GRanges object
structural_variants_genes_samples <- as.data.frame(sv_annotated) %>%
  dplyr::select(gene_name, Sample) %>%
  dplyr::group_by(gene_name, Sample) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()