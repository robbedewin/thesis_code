# Generates an HTML report summarizing the somatic variant calls using the maftools R package.
rule maftools_analysis:
    input:
        maf="results/mutect2/P013/P013_annotated_lifted_variants.maf",
    output:
        html="results/maftools/maftools_summary.html",
        pdf="results/maftools/maftools_oncoplot.pdf",
    params:
        dir="results/mutect2/"
    log:
        "logs/maftools/maftools_analysis.log"
    script:
        "scripts/maftools.R"
    