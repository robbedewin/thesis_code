.libPaths()
library(maftools)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(data.table)
library(tidyverse)
library(R.utils)
library(BSgenome.Hsapiens.UCSC.hg38)
library(grid)
library(ggrepel)
library(stringr)
library(rtracklayer)
.libPaths(c("/home/rstudio/host/lig/home/rcools/Rlibs/rocker-tidyverse/4.2.2", .libPaths()))
library(StructuralVariantAnnotation)

##Function to list all the files we need for the analysis:

list_files <- function(OUTBASE, pattern, samples_to_exclude,pattern_files) {
  dir_list <- list.files(path = OUTBASE, pattern = "PTCL", full.names = FALSE, recursive = FALSE)
  #Print the list of directries
  print(dir_list)
  #exclude samples if you want
  if (!is.null(samples_to_exclude) && length(samples_to_exclude) > 0) {
    dir_list <- dir_list[!dir_list %in% samples_to_exclude]
  }
  print(length(dir_list))
  #Create the full paths for each directory
  paths <- paste0(OUTBASE, "/", dir_list)
  print(paths)
  #List the files matching a specific pattern
  # List the files matching the pattern_files in the specified paths
  maffiles <- list.files(path = paths, pattern = pattern_files, full.names = TRUE, recursive = TRUE)
  # Print the list of matched files
  print(maffiles)
  # Print the length of the list of matched files
  print(length(maffiles))
  
  # Return the list of matched files
  return(maffiles)
  
}


samples_to_exclude <- c("PTCL31_T", "PTCL27_T", "PTCL20_T")
###List the maf files (SNV)
OUTBASE <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/SNV_calling/Mutect2/new_PON"
pattern <- "PTCL"
pattern_files <- ".maf$"
maffiles <- list_files(OUTBASE, pattern, pattern_files = pattern_files, samples_to_exclude = samples_to_exclude)

###List the ascat files (CNA)
OUTBASE <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/CNA_calling/ASCAT_data/matched_mode"
pattern <- "PTCL"
pattern_files <- ".segments.txt$"
samples_to_exclude <- c("PTCL31", "PTCL27", "PTCL20")
segfiles <- list_files(OUTBASE, pattern, pattern_files = pattern_files, samples_to_exclude = samples_to_exclude)

##List the GRIDSS files (SV)
OUTBASE <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/SV_calling/GRIDSS/calls/sv/matched_mode"
pattern <- "PTCL"
pattern_files <- "T.somatic.vcf.bgz$"
samples_to_exclude <- c("PTCL31_T", "PTCL27_T", "PTCL20_T")
SVABA_files <- list_files(OUTBASE, pattern, pattern_files = pattern_files, samples_to_exclude = samples_to_exclude)

##Start the SNV plot to get an overview:
y <- plotmafSummary(maf_files_matched_cohort, showBarcodes = TRUE)
x <- oncoplot(maf = maf_files_matched_cohort, top = 20, fontSize = 0.4, showTumorSampleBarcodes = T, writeMatrix = F)

##Make a dataframe and keep the columns which are relevant for the downstream analysis
dataframe_SNV <- maf_files_matched_cohort@data
dataframe_SNV <- as.data.frame(dataframe_SNV)
dataframe_SNV <- dataframe_SNV %>% dplyr::select(Hugo_Symbol, Chromosome, Start_Position,
                                          End_Position, Variant_Classification, Variant_Type,
                                          Tumor_Sample_Barcode, tumor_f,Protein_Change)

##Start looking at the ASCAT data, keep the regions where there is high-level amplification or LOH
#Load in the QC-data as-well to alter this analysis for the samples who are WGD:
ASCAT_data <- lapply(segfiles, function(x){
  y <- fread(x)
})

ASCAT_data <- lapply(ASCAT_data, function(x) {
  x$chr <- gsub("^([0-9, X]+)$", "chr\\1", as.character(x$chr))
  return(x)  # Return the modified data table
})

# Read in QC-file
qc <- file.path("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/CNA_calling/ASCAT_data/QC_data_entire_cohort.csv")
qc <- fread(qc)
qc_WGS <- subset(qc, select = c(Sample, QC.purity, QC.WGD))
WGD_subset <- qc_WGS[qc_WGS$QC.WGD ==1,]

### Make a function to spot LOH(losses), AMP and HOM_DEL

add_mutation_type_CN <- function(x, WGD_subset) {
  # Extract the sample name from x
  sample_name <- unique(x$sample)
  
  # Check if the sample is in WGD_subset
  if (!(sample_name %in% WGD_subset$Sample)) {
    x$mutation_type <- ifelse(x$nMinor < 1, "losses",
                              ifelse(x$nMinor == 0 & x$nMajor == 0, "Gene_homo_deletion",
                                     ifelse(x$nMajor >= 4, "Gene_amplification", NA)
                              ))
  } else {
    x$mutation_type <- ifelse(x$nMinor < 1, "losses",
                              ifelse(x$nMinor == 0 & x$nMajor == 0, "Gene_homo_deletion",
                                     ifelse(x$nMajor >= 6, "Gene_amplification", NA)))
  }
  return(x)
}


ASCAT_data_annotated <- lapply(ASCAT_data, function(x){
  y <- add_mutation_type_CN(x,WGD_subset)
})

###Annotate the regions with gene names: 
GFFfile <- readGFF("/home/rstudio/host/lig/references/GRCh38.alt-masked-V2-noALT/annotation/gencode.v43.basic.annotation.gtf")
GFFfile <- GFFfile[GFFfile$type == "gene", ]
GFFfile <- GFFfile[ , c(1,4,5,11)]
setDT(GFFfile)
##Remove Na's first
ASCAT_data_annotated_only <- lapply(ASCAT_data_annotated, function(x){
  na.omit(x)
})

#Find Overlap with the data.table
lapply(ASCAT_data_annotated_only,function(x){
  setDT(x)
})

setkey(GFFfile, seqid, start, end)
lapply(ASCAT_data_annotated_only, function(x){
  setkey(x,chr, startpos,endpos)
})

# Get overlaps. Make sure that start and end does not contain NA values and are integers/numeric
overlaps = lapply(ASCAT_data_annotated_only, function(x){
  foverlaps(x, GFFfile)
})

#Read in the COSMIC dataframe and Adjust in the format needed for the analysis:
cosmic <- fread("/home/rstudio/host/lig/references/COSMIC_CANCER_GENES_v98.tsv")
names(cosmic)
head(cosmic,5)
cosmic_filtered <- cosmic %>% dplyr::select('Gene Symbol','Role in Cancer', 'Genome Location')

COSMIC_data_frame <- cosmic_filtered %>%
  separate('Genome Location', into = c("Chromosome", "Position"), sep = "\\:") %>%
  separate(Position, into = c("Start_Position", "End_Position"), sep = "\\-")

COSMIC_data_frame$Chromosome <- paste0("chr",COSMIC_data_frame$Chromosome)
#First do a general overlap with data.table and filter later:
ALL_ASCAT <- do.call(rbind,overlaps)
ALL_ASCAT <- ALL_ASCAT %>% dplyr::select(chr,start,end,gene_name,sample,mutation_type)
setDT(ALL_ASCAT)
setDT(COSMIC_data_frame)

COSMIC_data_frame <- COSMIC_data_frame %>% rename(gene_name = 'Gene Symbol')

ALL_ASCAT <- na.omit(ALL_ASCAT)
COSMIC_data_frame <- na.omit(COSMIC_data_frame)
COSMIC_data_frame[, Start_Position := as.integer(Start_Position)]
COSMIC_data_frame[, End_Position := as.integer(End_Position)]

setkey(ALL_ASCAT, gene_name, chr, start, end)
setkey(COSMIC_data_frame, gene_name,Chromosome, Start_Position, End_Position)

COSMIC_ASCAT <- foverlaps(ALL_ASCAT, COSMIC_data_frame)
##Clean up the dataframe:
keep_COSMIC_ASCAT <- na.omit(COSMIC_ASCAT)
##Keep the TSG in losses and the oncogenes in the gained regions
keep_COSMIC_ASCAT <- keep_COSMIC_ASCAT %>%
  dplyr::select(gene_name, chr, 'Role in Cancer', Start_Position, End_Position, mutation_type, sample) %>%
  rename(Role = 'Role in Cancer') %>%
  dplyr::filter((mutation_type == "Gene_amplification" & str_detect(Role, "oncogene")) |
                  (mutation_type == "losses" & str_detect(Role, "TSG")))

##Additional filtering: keep TSG if there is a coding SNV in that gene
setDT(dataframe_SNV)
setDT(keep_COSMIC_ASCAT)

setkey(dataframe_SNV, Tumor_Sample_Barcode,Hugo_Symbol,Chromosome,Start_Position,End_Position)
setkey(keep_COSMIC_ASCAT, sample, gene_name,chr,Start_Position,End_Position)

SNV_COSMIC_ASCAT <- foverlaps(keep_COSMIC_ASCAT, dataframe_SNV)
SNV_COSMIC_ASCAT <- na.omit(SNV_COSMIC_ASCAT)
##All overlaps are in LOH regions -> 8 coding high impact TSG
oncogenes <- keep_COSMIC_ASCAT %>% dplyr::filter(mutation_type == "Gene_amplification")
names(oncogenes)
TSV <- SNV_COSMIC_ASCAT %>% dplyr::select(gene_name,chr,Role,Start_Position,End_Position,mutation_type,sample)

oncogens_and_TSV <- rbind(oncogenes,TSV)

###DONE##

#Let's do the last layer:SVs -> For this analysis use the GRIDSS calls, widen and filter:
sample_name <- basename(dirname(SVABA_files))

SVABA_dataframe <- lapply(X=SVABA_files, FUN = function(x) {
    vcfs <- VariantAnnotation::readVcf(x, genome = 'hg38')
    gr <- c(breakpointRanges(vcfs),breakendRanges(vcfs))
    wide <- resize(gr, width = 20000, fix = "center")
    return(wide)
  })

SVABA_dataframe <- lapply(SVABA_dataframe, function(x){
 y <- as.data.frame(x)
})

for (i in seq_along(SVABA_dataframe)) {
  SVABA_dataframe[[i]]$Tumor_Sample_Barcode <- sample_name[i]
}

SVABA_dataframe_all <- do.call(rbind,SVABA_dataframe)
#Filter and keep the columns we need
SVABA_filtered <- SVABA_dataframe_all %>% dplyr::select(seqnames,start,end,QUAL,svtype,svLen,Tumor_Sample_Barcode) %>%
  dplyr::filter(QUAL >= 1000)

SVABA_filtered <- na.omit(SVABA_filtered )
##Annotate with gene names:
setkey(GFFfile, seqid, start, end)
setDT(SVABA_filtered)
setkey(SVABA_filtered , seqnames,start,end)

SVABA_annotated <- foverlaps(SVABA_filtered,GFFfile)
SVABA_fill_and_annot  <- na.omit(SVABA_annotated)
##Done##
#Combine the data and first look at the cancer (COSMIC) associated genes:
oncogens_and_TSV 
COSMIC_genes <- COSMIC_data_frame$gene_name
COSMIC_SVs <- SVABA_fill_and_annot %>% dplyr::filter(gene_name %in% COSMIC_genes)
COSMIC_SNVs <- dataframe_SNV %>% dplyr::filter(Hugo_Symbol %in% COSMIC_genes)

##Combine all the impactfull mutations in one dataframe and start the plotting:
COSMIC_impact_SNVS 
oncogens_and_TSV 
COSMIC_SVs
##Rename so we have the same columns and column names:
oncogens_and_TSV_rename <- oncogens_and_TSV %>% dplyr::select(-Role) %>%
  rename(start=Start_Position,end=End_Position)


COSMIC_impact_SNVS_rename <- COSMIC_impact_SNVS %>% dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode) %>%
  rename(gene_name = Hugo_Symbol, chr = Chromosome, mutation_type = Variant_Classification, sample = Tumor_Sample_Barcode, start=Start_Position,end=End_Position)
  
SVABA_fill_and_annot_rename <- COSMIC_SVs %>% dplyr::select(seqnames, start,end,gene_name,svtype,Tumor_Sample_Barcode) %>%
  rename(chr = seqnames, mutation_type = svtype, sample = Tumor_Sample_Barcode)


ALL_COSMIC_IMPACT_MUT <- rbind(oncogens_and_TSV_rename ,COSMIC_impact_SNVS_rename,SVABA_fill_and_annot_rename,CH_filtered)

getwd()

##See which genes are to mutated in the dataframe:
all <- as.data.frame(with(ALL_COSMIC_IMPACT_MUT, table(gene_name)))
all_ordered <- all %>% arrange(desc(Freq))
top <- head(all_ordered, 50)

##First do the lotting for the COSMIC genes and then Add the clonal Hemato genes:
ALL_COSMIC_IMPACT_MUT <- ALL_COSMIC_IMPACT_MUT %>% mutate(mutation_type = recode(mutation_type, 
                                                                                 "In_Frame_Del" = "In_Frame_mut",
                                                                                 "In_Frame_Ins" = "In_Frame_mut"))
unique(ALL_COSMIC_IMPACT_MUT$mutation_type)
ALL_COSMIC_IMPACT_MUT <- ALL_COSMIC_IMPACT_MUT %>% mutate(mutation_type = recode(mutation_type, 
                                                                                 "Frame_Shift_Del" = "Frame_Shift_mut",
                                                                                 "Frame_Shift_Ins" = "Frame_Shift_mut"))

for_use <- ALL_COSMIC_IMPACT_MUT %>% dplyr::select(sample,mutation_type,gene_name)
for_use <- for_use %>% dplyr::filter(gene_name %in% top$gene_name)
unique(for_use$mutation_type)
wide_df <- pivot_wider(for_use, names_from = sample, values_from = mutation_type, values_fn = list(mutation_type = toString))
wide_df[is.na(wide_df)] <- ""
wide_df
wide_df <- column_to_rownames(wide_df, var = "gene_name")
wide_df <- wide_df %>%
  mutate_all(~ str_replace_all(., ",", ";"))

##Start the plotting with oncorint from Complex Heatmap:
library(RColorBrewer)
display.brewer.all(colorblindFriendly = F)
brewer.pal(n = 24, name = "Set3")

#Make Heatmap
col = c("losses" = "#4E84C4", "Gene_amplification" = "#D16103", "HOM_del" = "#F4EDCA", 
        "BND" = "#FFDB6D", "Nonsense_Mutation" = "#52854C", 
        "Missense_Mutation" = "#FC8D62", "Frame_Shift_mut" = "#293352", 
        "Splice_Site" = "#999999", "In_Frame_mut"= "#66C2A5", "Frame_Shift_Ins"= "#B3B3B3", 
        "In_Frame_Ins"="#66C2A5", "Nonstop_Mutation"=  "#FFD92F", "Precursor/CH-SNV" = "#BC80BD"
        )
# Create functions for colors present in 'col'
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  losses = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), 
              gp = gpar(fill = col["losses"], col = NA))
  },
  Gene_amplification = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), 
              gp = gpar(fill = col["Gene_amplification"], col = NA))
  },
  HOM_del = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), 
              gp = gpar(fill = col["HOM_del"], col = NA))
  },
  BND = function(x, y, w, h) {
    grid.polygon(unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
                 unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
                 gp = gpar(fill = col["BND"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x - 0.25*w, y, 0.5*w - unit(0.5, "pt"), h - unit(0.5, "pt"),
                 gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x + 0.25*w, y, 0.5*w - unit(0.5, "pt"), h - unit(0.5, "pt"),
                 gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Frame_Shift_mut = function(x, y, w, h) {
    grid.rect(x, y + 0.25*h, w - unit(0.5, "pt"), 0.5*h - unit(0.5, "pt"),
                 gp = gpar(fill = col["Frame_Shift_mut"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
      grid.rect(x, y - 0.25*h, w - unit(0.5, "pt"), 0.5*h - unit(0.5, "pt"), 
                 gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  In_Frame_mut = function(x, y, w, h) {
    grid.rect(x, y - 0.25*h, w - unit(0.5, "pt"), 0.5*h - unit(0.5, "pt"), 
              gp = gpar(fill = col["In_Frame_mut"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y + 0.25*h, w - unit(0.5, "pt"), 0.5*h - unit(0.5, "pt"),
                 gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y - 0.25*h, w - unit(0.5, "pt"), 0.5*h - unit(0.5, "pt"),
                 gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  Nonstop_Mutation  = function(x, y, w, h) {
    grid.polygon(unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
                 unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
                 gp = gpar(fill = col["Nonstop_Mutation "], col = NA))
  },
  'Precursor/CH-SNV' = function(x, y, w, h) {
    grid.polygon(unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
                 unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
                 gp = gpar(fill = col["Precursor/CH-SNV"], col = NA))
  }
)
library(ComplexHeatmap)


column_title = "Mutational landscape of PTCL-NOS COSMIC genes"
heatmap_legend_param = list(title = "Type of mutations", at = c("losses", "Gene_amplification", "HOM_del", "BND","Nonsense_Mutation", "Missense_Mutation", "Frame_Shift_mut", "Splice_Site", "In_Frame_mut", "Nonstop_Mutation",'Precursor/CH-SNV'),
                            labels = c("LOH", "High-level amp", "HOM_del", "BND", "Nonsense_Mutation", "Missense_Mutation", "Frame_Shift_mut", "Splice_Site", "In_Frame_mut", "Nonstop_Mutation",'Precursor/CH-SNV'))
p1 <- oncoPrint(wide_df ,
                alter_fun = alter_fun,show_column_names = T,col = col,alter_fun_is_vectorized = FALSE,
                column_title = column_title, heatmap_legend_param = heatmap_legend_param,pct_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 7), remove_empty_rows = F,remove_empty_columns = F, row_names_gp = gpar(fontsize = 5))
p1q