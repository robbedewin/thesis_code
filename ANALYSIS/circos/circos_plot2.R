install.packages("circlize")
library(circlize)


library(data.table)
library(dplyr)

# Load MAF, VCF, and Segment data
maf <- fread("/staging/leuven/stg_00096/home/rdewin/WGS/results/vcf2maf/P011/P011_mutect_pass_variants_annotated.maf")
segments <- fread("/staging/leuven/stg_00096/home/rdewin/WGS/results/ascat/P011/P011_tumor.segments.txt")
vcf <- fread("/staging/leuven/stg_00096/home/rdewin/WGS/results/gridss/P011/P011_high_confidence_somatic.vcf",
             skip = "#CHROM",  # Start reading from the header line
             fill = TRUE)      # Handle missing columns in some rows
   



maf_data <- maf %>%
  dplyr::select(Chromosome = Chromosome, Start = Start_Position, End = End_Position, Variant = Variant_Classification)


vcf_data <- vcf %>%
  filter(grepl("BND|DEL|DUP|INV", INFO)) %>%
  mutate(Chromosome = gsub("^chr", "", "#CHROM"), 
         Start = POS, 
         End = as.numeric(gsub(".*END=([0-9]+).*", "\\1", INFO))) %>%
  dplyr::select(Chromosome, Start, End, SVTYPE = INFO)


segments_data <- segments %>%
  mutate(Chromosome = gsub("^chr", "", chr)) %>%
  dplyr::select(Chromosome, Start = startpos, End = endpos, nMajor, nMinor)

chrom_info <- read.table("/staging/leuven/stg_00096/home/rdewin/ANALYSIS/circos/chrom_sizes.txt", 
                         header = TRUE, 
                         col.names = c("chr", "start", "end"))


circos.genomicInitialize(chrom_info)

circos.clear()
circos.initializeWithIdeogram(species = "hs1")

circos.genomicTrack(maf_data, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = "red", cex = 0.5, pch = 16)
})

circos.genomicTrack(vcf_data, panel.fun = function(region, value, ...) {
  circos.genomicLinks(region, value, col = "blue", lwd = 2)
})

circos.genomicTrack(segments_data, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = ifelse(value$nMajor > 2, "green", "orange"), border = NA)
})

circos.clear()
circos.initializeWithIdeogram(species = "chm13", plotType = c("axis", "labels"))


