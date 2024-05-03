library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(GenomicRanges)
library(RColorBrewer)
library(parallel)
library(grid)
library(ggrepel)
library(stringr)
library(tidyverse)

# Open the log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#Define all the inutfiles
sample_id <- snakemake@wildcards[["sample"]]

outbase <- snakemake@output[[1]]


#List all files for ASCAT (some files are tumour only, others not)
OUTBASE <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data"
sampleid <- c("PTCL1","PTCL23", "PTCL15", "PTCL17", "PTCL23", "PTCL26", "PTCL27", "PTCL3", "PTCL5", "PTCL9", "PTCL11","PTCL14", "PTCL18", "PTCL28")
OUTBASE1 <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/Tumour_only"
sampleid_tumour_only <- c("PTCL12", "PTCL13", "PTCL16")
OUTBASE2 <- "/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/Tumour_only/new_update"
sampleid_tumour_only_new_update <- c("PTCL19", "PTCL20", "PTCL21", "PTCL22", "PTCL25")

getwd()
setwd("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/")
?list.dirs
#Normal
dir_list <- list.dirs(path = OUTBASE, full.names = FALSE)
PTCL_dirs <- dir_list[dir_list %in% sampleid]
PTCL_dirs
paths <- paste0(OUTBASE, "/",PTCL_dirs)
paths
segfiles <- list.files(path = paths, pattern = ".segments.txt$", full.names = T, recursive = T)
segfiles
#Tumour_only
dir_list1 <- list.dirs(path = OUTBASE1, full.names = FALSE)
dir_list1
PTCL_dirs1 <- dir_list1[dir_list1%in% sampleid_tumour_only]
PTCL_dirs1
paths1 <- paste0(OUTBASE1, "/",PTCL_dirs1)
paths1
segfiles1 <- list.files(path = paths1, pattern = ".segments.txt$", full.names = T, recursive = T)
segfiles1
#Tumour_only_update
dir_list2 <- list.dirs(path = OUTBASE2, full.names = FALSE)
dir_list2
PTCL_dirs2 <- dir_list2[dir_list2%in% sampleid_tumour_only_new_update]
PTCL_dirs2
paths2 <- paste0(OUTBASE2, "/",PTCL_dirs2)
paths2
segfiles2 <- list.files(path = paths2, pattern = ".segments.txt$", full.names = T, recursive = T)
segfiles2
#combine all seg files
all_ascat_seg_files <- c(segfiles,segfiles1,segfiles2)
all_ascat_seg_files

#list QC files
qcfiles <- list.files(path = paths, pattern = "objects.Rdata$", full.names = T, recursive = T)
qcfiles
qcfiles1 <- list.files(path = paths1, pattern = "_objects.Rdata$", full.names = T, recursive = T)
qcfiles1
qcfiles2 <- list.files(path = paths2, pattern = "_objects.Rdata$", full.names = T, recursive = T)
qcfiles2
all_qc_files23 <- c(qcfiles,qcfiles1,qcfiles2)
all_qc_files23

# Initialize a data frame to store the "QC" components
qc_data <- data.frame()

# Loop through each file
for (x in all_qc_files23) {
  # Load the object from the file
  load(x)
  
  # Extract the "QC"
  qc_component <- QC  
  
  # Create a data frame row with the "QC" component and the file name
  row <- data.frame(Sample = x, QC = qc_component)
  
  # Add the row to the qc_data data frame
  qc_data <- rbind(qc_data, row)
}
# Write the qc_data data frame to a tab-separated .csv file
write.table(qc_data, "qc_data.csv", sep = "\t", row.names = FALSE, quote = FALSE)
qc_file <- paste("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/qc_data.csv")
getwd()
qcdata <- read.delim("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/qc_data.csv", as.is = T)
qcdata

#Read in segfiles and create a GR with "annotations"
segdata <- lapply(X = all_ascat_seg_files, FUN = function(x) {
  y <- read.delim(x, as.is = T)
  y$chr <- gsub ("^([0-9, X]+)$", "chr\\1", as.character(y$chr))
  y <- GRanges(seqnames = y$chr, ranges = IRanges(start = y$startpos, end = y$endpos),
               sample = y$sample, nMajor = y$nMajor, nMinor = y$nMinor, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
  y <- keepSeqlevels(x = y, value = paste0("chr", 1:22), pruning.mode = "coarse")
  return(y)
})
segdata[[1]]
length(segdata)
names(segdata) <- sapply(X = segdata, FUN = function(x) x$sample[1])
#qcdata$Sample <- sapply(X = segdata, FUN = function(x) x$sample[1])
qcdata

#Overal QC
qcdata$
#Plotting loh and ploidy
m <- (2.9 - 1)/(0-0.93)
p1 <- ggplot(data = qcdata, mapping = aes(x = QC.LOH, y = QC.ploidy, colour = QC.WGD))
p1 <- p1 + geom_point() + geom_abline(slope = m, intercept = 2.9) + theme_minimal() + theme(plot.background = element_rect(fill = "white"))
p1
setwd("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/plots/")
ggsave(filename = "frac_loh_vs_ploidy.png", plot = p1, width = 10, height = 7 )

p1 <- ggplot(data = qcdata, mapping = aes(x = QC.LOH, y = QC.ploidy, colour = QC.WGD))
p1 <- p1 + geom_point() + geom_abline(slope = m, intercept = 2.9) + theme_minimal() + theme(plot.background = element_rect(fill = "white"))
p1 <- p1 + geom_label_repel(aes(label = Sample), label.size  = 0.3, size = 1.8, max.overlaps = Inf)
p1
ggsave(filename = "annotated_frac_loh_vs_ploidy.png", plot = p1, width = 10, height = 7 )

#frac genome 
checkfiles <- sapply(X = segfiles, FUN = function(x) file.size(x) > 1)
segfiles <- segfiles[checkfiles]
fracaberdip <- sapply(X = segdata, FUN = function(x) sum(width(x[which(x$nMajor != 1 & x$nMinor != 1), ])) / sum(width(x)))
fracabertet <- sapply(X = segdata, FUN = function(x) sum(width(x[which(x$nMajor != 2 & x$nMinor != 2), ])) / sum(width(x)))
fracloh <- sapply(X = segdata, FUN = function(x) sum(width(x[which(xor(x = x$nMinor == 0, y = x$nMajor == 0)), ])) / sum(width(x)))
outdf1 <- data.frame(sampleid = sub(pattern = ".segments.txt", replacement = "", x = basename(all_ascat_seg_files)),
                     frac_aber_dip = fracaberdip, frac_aber_tet = fracabertet, frac_loh = fracloh)

outdf1
qcdata$QC.GI
qcdata
#Use of pipe statement instead of the or statement

sum(width(segdata$PTCL11_T[which(segdata$PTCL11_T$nMajor != 1 & segdata$PTCL11_T$nMinor != 1), ])) / sum(width(segdata$PTCL11_T))
sum(width(segdata$PTCL11_T[which(segdata$PTCL11_T$nMajor != 1 | segdata$PTCL11_T$nMinor != 1), ])) / sum(width(segdata$PTCL11_T))
qcdata
#add values to qcdata, not the same as GI column...
qcdata$QC.frac_aber <- ifelse(qcdata$QC.WGD, outdf1$frac_aber_tet, outdf1$frac_aber_dip)
qcdata$QC.frac_aber[is.na(qcdata$ploidy)] <- NA
qcdata
MYC_pos <- c("PTCL12_T", "PTCL16_T", "PTCL19_T", "PTCL22_T", "PTCL25_T", "PTCL26_T", "PTCL28_T", "PTCL5_T")
qcdata$MYC_sig <- ifelse(qcdata$Sample %in% MYC_pos,1,0)
qcdata

colpal <- colorRampPalette(brewer.pal(12, "Set2"))
myPal <- colpal(length(unique(qcdata$Sample)))

#Plotting
p1 <- ggplot(data = qcdata, aes(x = Sample, y = QC.frac_aber, size = 0.2, color = MYC_sig)) + 
  geom_jitter(width = 0.25, show.legend = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p1
ggsave(filename = "Frac_aberatted_with_Myc.png", plot =p1, width=10, height = 5)

p1 <- ggplot(data = qcdata, aes(x = Sample, y = QC.GI, size = 0.2, color = MYC_sig)) + 
  geom_jitter(width = 0.25, show.legend = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p1
ggsave(filename = "QC.GI_of_the_cohort.png", plot = p1, width = 10, height = 5)


p2 <- ggplot(data = qcdata, mapping = aes(x = Sample, y = QC.GI, size = 0.2)) + 
  geom_jitter(mapping = aes(color = Sample), width = 0.25, show.legend = TRUE)
p2 <- p2 + theme(axis.text.x = element_text(size = 12))  # Adjust the size (e.g., size = 12)
p2 <- p2 + scale_fill_manual(values = myPal) + scale_color_manual(values = myPal) + theme_minimal() + theme(plot.background = element_rect(fill = "white")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p2
ggsave(filename = "QC.GI.png", plot = p2, width = 10, height = 5)

segdata
#Make overal ASCAT plot
length(segdata)
segdata
sampleseg <- unlist(GRangesList(segdata[qcdata$Sample]))
sampleseg
joints <- disjoin(sampleseg)
joints
list(bins = joints)
out <- mapply(sampleid = qcdata$Sample, wgd = qcdata$QC.WGD, nsamples = 1, MoreArgs = list(bins = joints),
              FUN = function(sampleid, wgd, nsamples, bins){
                if (wgd == 0) {
                  gains <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMajor > 1])/nsamples
                  losses <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMinor < 1])/nsamples
                } else {
                  gains <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMajor > 2])/nsamples
                  losses <- countOverlaps(query = bins, subject = sampleseg[sampleseg$sample == sampleid & sampleseg$nMinor < 2])/nsamples
                }
                return(list(gains=gains, losses =losses))
              }, SIMPLIFY = F)

joints$ngains <- Reduce(f = '+', x = lapply(X = out, FUN = function(x) x$gains))
joints$nlosses <- Reduce(f = '+', x = lapply(X = out, FUN = function(x) x$losses))
joints$ncases <- 20

joints$gains <- joints$ngains / joints$ncases
joints$losses <- joints$nlosses / joints$ncases

cumdist <- setNames(object = c(0, cumsum(as.numeric(seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:22]))), nm = paste0("chr", 1:23))
joints$cumstart <- start(joints) + cumdist[as.character(seqnames(joints))]
joints$cumend <- end(joints) + cumdist[as.character(seqnames(joints))]
joints

p1 <- ggplot(data = as.data.frame(joints)) + geom_rect(mapping = aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = gains), fill = "#fc8d59", alpha = .6)
p1 <- p1 + geom_rect(mapping = aes(xmin = cumstart, xmax = cumend, ymin = 0, ymax = -losses), fill = "#91bfdb", alpha = .6)
p1 <- p1 + geom_vline(xintercept = cumdist) + ylim(c(-1,1))
p1 <- p1 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
p1 <- p1 + scale_x_continuous(breaks = filter(x = cumdist, filter = c(.5, .5), sides = 1)[-1], labels = 1:22) + theme_minimal()+ theme(plot.background = element_rect(fill = "white"))
p1 <- p1 + labs(x = "", y = "Frequency")
p1
ggsave(filename = "ascat_overview.png", plot = p1, width =10, height= 5)


#Additional QC from the ASCAT plots
#QC was determined as follows:
#  - **No ASCAT fit**: ASCAT could not find an optimal ploidy and purity value.
#- **Likely normal**: CNA profile matches with a normal sample. Filter: LOH<0.1 & GI<0.1 & purity=1 & ((sex='XX' & ploidy>=1.99 & ploidy<=2.01) | (sex='XY' & ploidy>=1.945 & ploidy<=1.965)).
#- **Noisy logR**: logR track (either tumour or normal) is noisy. Filter: tumour_mapd>=0.75 | normal_mapd>=0.75.
#- **Likely wrong germline**: matched normal is likely to come from another patient. Filter: tumour_mapd>0.4 & normal_mapd>0.4 & frac_homo>0.1.
#- **Oversegmented**: difference in segmentation between logR and BAF tracks. Filter: n_segs_logRBAF_diff>250.
#- **Wrong fit**: CNA profile is incorrect, with extreme losses. Filter: mode_majA=0.
#- **HD size**: CNA profile contains large segments with homozygous deletion. Filter: homdel_largest>=20e6 | homdel_size>=40e6.
#- **Contamination or swap** and **Large CNV**: germline data contains copy-number changes. Filter: we ran aspcf on all germline samples and manually reviewed cases where events were spotted. Cases with large stretches of homozygosity were considered as valid and were not discarded.
#- **Pass**: sample has correct metrics.

qcdata$Likely_normal <- ifelse(qcdata$QC.LOH < 0.1 & qcdata$QC.GI < 0.1 & qcdata$QC.purity == 1 & ((qcdata$QC.sex == 'XX' & qcdata$QC.ploidy >= 1.99 & qcdata$QC.ploidy <= 2.01) | (qcdata$QC.sex == 'XY' & qcdata$QC.ploidy >= 1.945 & qcdata$QC.ploidy <= 1.965)), TRUE, FALSE)
qcdata                               
qcdata$Noisy_LogR <- ifelse(qcdata$QC.tumour_mapd >= 0.75 | (!is.na(qcdata$QC.normal_mapd) & qcdata$QC.normal_mapd >= 0.75), TRUE, FALSE)
qcdata$Likely_wrong_Germline <- ifelse(qcdata$QC.tumour_mapd > 0.4 | (!is.na(qcdata$QC.normal_mapd) & qcdata$QC.normal_mapd >= 0.4) & qcdata$QC.frac_homo > 0.1, TRUE, FALSE)
qcdata$Oversegmented <- ifelse(qcdata$QC.n_segs_logR-qcdata$QC.n_segs_BAF > 250, TRUE, FALSE)
qcdata$wrong_fit <- ifelse(qcdata$QC.mode_majA == 0, TRUE, FALSE)
qcdata$HD_size <- ifelse(qcdata$QC.homdel_largest >= 20e6 | qcdata$QC.homdel_size >= 40e6, TRUE, FALSE)
write.table(qcdata, "qc_data.csv", sep = "\t", row.names = FALSE, quote = FALSE)
qcdata

p1 <- ggplot(data = qcdata, mapping = aes(x = QC.n_segs_BAF, y = QC.n_segs_logR, colour = QC.WGD))
p1 <- p1 + geom_point() + geom_abline(slope = m, intercept = 250) + theme_minimal() + theme(plot.background = element_rect(fill = "white"))
p1 <- p1 + geom_label_repel(aes(label = Sample), label.size  = 0.3, size = 1.8, max.overlaps = Inf)
p1
ggsave(filename = "oversegmented.pdf",plot = p1, width =10, height= 5 )
setwd("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/plots")

m <- 1
m <- (2.9 - 1)/(0-0.93)

# Making Heatmap of the CNV landscape (following some steps of cnSpec)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
library(circlize)

#GRanges of all chromosomes and divide in in windows of 1MB
chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_gr
chr_window = makeWindows(chr_gr, w = 5e6)
chr_window
#joints$cumstart <- start(joints) + cumdist[as.character(seqnames(joints))]
#joints$cumend <- end(joints) + cumdist[as.character(seqnames(joints))]

new_window <- GRanges(seqnames = seqnames(joints), ranges = ranges(joints))

segdata[[1]]
#calculate the average signals in the 1MB windows 
average_in_window = function(window, gr, v, empty_v=NA, method = "weighted") { 
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    }
  }
  
  return(u)
}
sampleseg$total_CN <- sampleseg$nMajor + sampleseg$nMinor
PTCL1_test <- sampleseg[sampleseg$sample == "PTCL1_T",]
#sampleseg
#CN_data <- average_in_window(new_window,PTCL1_test,PTCL1_test$total_CN)

CNsegdata <- lapply(segdata, FUN = function(x){
  x$total_CN <- x$nMinor + x$nMajor
  return(x)
})
CNsegdata


segdata
total_CN_data <- lapply(X = CNsegdata, FUN = function(x) {
  y <- average_in_window(new_window, x,x$total_CN)
  return(y)
})
total_CN_data[[20]]
names(total_CN_data) <- sapply(X = segdata, FUN = function(x) x$sample[1])

chr_window_df <- as.data.frame(new_window)
sample_names <- names(total_CN_data)

# Loop through each sample and add it as a new column
for (sample_name in sample_names) {
  chr_window_df[[sample_name]] <- total_CN_data[[sample_name]]
}
chr_window_df

segdata[[4]]

final_table <- chr_window_df %>% select(-width,-strand)
final_table

long_table <- final_table %>% pivot_longer(cols = starts_with("PTCL"), names_to="sample", values_to="cn")
long_table <- long_table %>% rename("chromosome" = "seqnames")
long_table <- long_table %>% rename("segmean" = "cn")
names(long_table)

#Complex Heatmap
col_fun = colorRamp2(c(0,2,6), c("blue","white","red"))
col_fun(seq(-3, 3))


final_table
joints


#GenVisR heatmap
BiocManager::install("GenVisR")
library(GenVisR)

genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg38",], max)
genomeBoundaries$chromStart <- 0
colnames(genomeBoundaries) <- c("chromosome", "end", "start")
cumdist

replace_na_with_2 <- function(x) {
  x[is.na(x)] <- 2
  return(x)
}
replace_na_with_2(long_table$segmean)
p1 <- cnSpec(long_table, y = genomeBoundaries, genome = "hg38", CNscale = "absolute")
p1
ggsave(filename = "heatmap_first_time.pdf",plot = p1, width =10, height= 5 )
setwd("/home/rstudio/host/lig/home/rcools/projects/data_Cools_2022/ASCAT_data/plots")
genomeBoundaries_2 <- genomeBoundaries
genomeBoundaries_2$start <- genomeBoundaries_2$start - 1e8
genomeBoundaries_2$end <- genomeBoundaries_2$end + 1e8
