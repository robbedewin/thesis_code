library("DESeq2")
library(scales)


# Open the log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")



parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

